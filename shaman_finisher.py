#!/usr/bin/env python
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html
from bioblend.galaxy import GalaxyInstance
import bioblend
import time
import os
import sys
import argparse
import json
import glob
import shutil
#import keyring
import zipfile
import smtplib
import socket
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        # If false then its a list
        if not isinstance(values, str):
            out = []
            for val in values:
                if os.path.isfile(val):
                    out += [os.path.abspath(os.path.expanduser(val))]
                elif os.path.isdir(val):
                    out += [os.path.abspath(os.path.expanduser(val)) + os.sep]
                else:
                    out += [val]
            setattr(namespace, self.dest, out)
        # Value is a string
        else:
            if os.path.isfile(values):
                setattr(namespace, self.dest,
                        os.path.abspath(os.path.expanduser(values)))
            elif os.path.isdir(values):
                setattr(namespace, self.dest,
                        os.path.abspath(os.path.expanduser(values)) + os.sep)

class galaxy:

    def __init__(self, task_file, done_dir, galaxy_url, galaxy_key, https_mode, message):
        #jobid
        #num_job

        self.galaxy_url = galaxy_url
        self.galaxy_key = galaxy_key
        self.gi = GalaxyInstance(url=galaxy_url, key=galaxy_key)
        self.gi.verify = https_mode
        self.task_file = task_file
        self.done_dir = done_dir
        self.message = message
        #self.num_job = num_job
        #self.jobid = jobid


    def load_json(self):
        """Load and validate Json
        """
        data_task_ok = None
        try:
            with open(self.task_file, "rt") as task:
                data_task = json.load(task)
                if isinstance(data_task, (list, tuple)):
                    data_task = data_task[0]
                #data_task["name"] = os.path.splitext(os.path.basename(self.task_file))[0]
                data_task_ok = data_task
        except IOError as err:
            err.extra_info("Error cannot open {0}".format(self.task_file)) 
        return data_task_ok


    def reconnect(self):
        """Reconnect to galaxy
        """
        self.gi = GalaxyInstance(url=self.galaxy_url, key=self.galaxy_key)
        self.gi.verify = False

    def zip_archive(self, list_downloaded_files, zip_file):
        """Extract tar file and build a clean zip file
        """
        #try:
        with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for file in list_downloaded_files:
                zipf.write(file, os.path.basename(file))
        #except IOError:
        #    self.logger.error("Error cannot open {0}".format(zip_file))

    def send_mail(self, result_file=None):
            """Send result by email
            """
            #try:
            fromaddr = "shaman@pasteur.fr"
            #bcc = ['amine.ghozlane@pasteur.fr']
            toaddr = self.data_task["mail"]
            msg = MIMEMultipart()
            msg['From'] = fromaddr
            msg['To'] = toaddr
            msg['Subject'] = "Shaman result"
            msg.attach(MIMEText(self.message, 'plain'))
            if result_file:
                if os.path.getsize(result_file) < 10000000:
                    part = MIMEBase('application', 'octet-stream')
                    with open(result_file, "rb") as attachment:
                        part.set_payload((attachment).read())
                        encoders.encode_base64(part)
                        part.add_header('Content-Disposition',
                                    "attachment; filename= {0}"
                                    .format(os.path.basename(result_file)))
                        msg.attach(part)
            if socket.gethostname() == "ShinyPro":
                server = smtplib.SMTP('smtp.pasteur.fr', 25)
                server.starttls()
                text = msg.as_string()
                toaddr = [toaddr] + ["amine.ghozlane@pasteur.fr"]
                server.sendmail(fromaddr, toaddr, text)
                server.quit()

    def download_result(self, history_id, list_result, result_dir):
        """Download tar archive from galaxy
        """
        list_downloaded_files = []
        success = True
        # Create output directory
        #try:
        if not os.path.isdir(result_dir):
            os.mkdir(result_dir)
        try:
            for result_type in list_result:
                for result_file in list_result[result_type]:
                    match = self.gi.histories.show_matching_datasets(
                                    history_id, result_file)
                    if len(match) >0:
                        res = result_dir + result_file + "." + result_type
                        self.gi.datasets.download_dataset(
                            match[0]['id'], file_path = res, 
                            use_default_filename=False, 
                            wait_for_completion=True, maxwait=60)
                        if os.stat(res).st_size == 0:
                            self.logger.error("File {0} is empty".format(res))
                            #success = False
                        else:
                            list_downloaded_files.append(res)
                    else:
                        self.logger.error("Match for result file: {} and result type: {} = {}".format(result_file, result_type, match))
            assert(len(list_downloaded_files) > 0)
        except bioblend.galaxy.datasets.DatasetTimeoutException:
            success = False
        except AssertionError:
            success = False
        return success, list_downloaded_files

    def run(self):
        """Upload, run galaxy workflow and dowload results
        """
        list_result = {}
        self.data_task = self.load_json()
        result_history = self.gi.histories.get_histories(name=self.data_task['result_history_name'])[0]
        print(result_history)
        # Load json data
        
        # Output
        zip_file = self.done_dir + os.sep + "shaman_" + self.data_task["name"].replace("file", "") + ".zip"
        result_dir = self.done_dir + os.sep + self.data_task["name"] + os.sep
        # Expected result files
        list_result['fasta'] = ["shaman_otu"]
        list_result['tsv'] = ["shaman_rdp_annotation", "shaman_otu_table",
                        "shaman_process_build", "shaman_process_annotation"]
        if self.data_task["type"] == "16S":
            list_result['biom'] = ["shaman_silva", "shaman_greengenes"]
        elif self.data_task["type"] == "18S" or self.data_task["type"] == "23S_28S":
            list_result['biom'] = ["shaman_silva"]
        else:
            list_result['biom'] = ["shaman_findley", "shaman_unite", "shaman_underhill"]
        # Add tree and annotation files
        list_result['tsv'] += [i + "_annotation"  for i in list_result['biom']]
        list_result['nhx'] = [i + "_tree"  for i in list_result['biom']]
        # Download results
        download_success, list_downloaded_files = self.download_result(
                                result_history['id'], list_result,
                                result_dir)
        #if os.path.isfile(result_file) and download_success:
        if download_success:
            # Prepare zip
            self.zip_archive(list_downloaded_files, zip_file)
            self.send_mail(zip_file)


def isdir(path):
    """Check if path is an existing file.
      Arguments:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file.".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def isfile(path):
    """Check if path is an existing file.
      Arguments:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def getArguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.add_argument('-u', dest='galaxy_url', type=str, #required=True,
                        default='https://galaxy.pasteur.fr',
                        #default='https://galaxy-dev.web.pasteur.fr',
                        #default='http://127.0.0.1:8080',
                        help='Url to galaxy (default https://galaxy.pasteur.fr).')
    parser.add_argument('-k', dest='galaxy_key', type=str, #required=True,
                        default='31f05d9edaa2228b66c538f43b0d5d52',
                        #default=keyring.get_password("galaxy", "aghozlan"),
                        #default='7ac30484f696937116f960531a05c2b6',
                        #default='f293dce7785a77c338db9e8b8df9922c',
                        help='User galaxy key (default 31f05d9edaa2228b66c538f43b0d5d52).')
    parser.add_argument('-i', dest='todo_file', type=isfile, required=True,
                        help='Todo job file.')
    #parser.add_argument('-j', dest='jobid', type=str, required=True,
    #                    help='Galaxy job id.')
    parser.add_argument('-m', dest='message', type=str, required=True,
                        help='Error message to send to user.')
    parser.add_argument('-w', dest='done_dir', type=isdir, required=True,
                        action=FullPaths, help='Path to the result directory.')
    parser.add_argument('-s', dest='https_mode', action='store_true',
                        default=False, help='Activate https verification.')
    args = parser.parse_args()
    return args
 
        


def main():
    """Main program
    """
    args = getArguments()
    djinn = galaxy(args.todo_file, args.done_dir, args.galaxy_url,
                   args.galaxy_key, args.https_mode, args.message)
    djinn.run()



if __name__ == '__main__':
    main()
