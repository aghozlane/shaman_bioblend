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
from threading import Thread
import daemon
import logging
from logging.handlers import RotatingFileHandler
import tempfile
import time
import os
import sys
import argparse
import json
import glob
import shutil
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders
import requests 



class galaxy(Thread):

    def __init__(self, logger, task_file, done_dir, error_dir, galaxy_url,
                 galaxy_key, num_job):
        Thread.__init__(self)
        self.logger = logger
        self.galaxy_url = galaxy_url
        self.galaxy_key = galaxy_key
        self.gi = GalaxyInstance(url=galaxy_url, key=galaxy_key)
        self.gi.verify = False
        self.task_file = task_file
        self.done_dir = done_dir
        self.error_dir = error_dir
        self.num_job = num_job
        self.dataset_ids = []

    def load_json(self):
        """Load and validate Json
        """
        data_task_ok = None
        try:
            try:
                with open(self.task_file, "rt") as task:
                    data_task = json.load(task)
                    if "paired" not in data_task:
                        raise ValueError("Paired information is missing in {0}".format(self.task_file))
                    isinstance(data_task["paired"], bool)
                    if data_task["paired"]:
                        for pair in ["path_R1", "path_R2"]:
                            if pair not in data_task:
                                raise ValueError(
                                    pair + " information is missing in {0}".format(self.task_file))
                            if not os.path.isdir(data_task[pair]):
                                raise ValueError(
                                    pair + " is not a file in {0}".format(self.task_file))
                        for pair in ["contaminant_R1", "contaminant_R2"]:
                            if pair not in data_task:
                                raise ValueError(
                                    pair + " information is missing in {0}".format(self.task_file))
                            if not os.path.isfile(data_task[pair]):
                                raise ValueError(
                                    pair + " is not a file in {0}".format(self.task_file))
                    else:
                        if "path" not in data_task:
                            raise ValueError(
                                "path information is missing in {0}".format(self.task_file))
                        if not os.path.isdir(data_task["path"]):
                            raise ValueError(
                                "path information is not a directory for {0}".format(self.task_file))
                        if "contaminant" not in data_task:
                            raise ValueError(
                                "path information is missing in {0}".format(self.task_file))
                        if not os.path.isfile(data_task["contaminant"]):
                            raise ValueError(
                                "contaminant information is not a directory for {0}".format(self.task_file))
                    if "host" not in data_task:
                        raise ValueError("host information is missing in {0}".format(self.task_file))
                    if "type" not in data_task:
                        raise ValueError("type information is missing in {0}".format(self.task_file))
                    if data_task["type"] not in ["16S_18S", "23S_28S",
                                                 "ITS", "WGS"]:
                        raise ValueError("type information is incorrect in {0}".format(self.task_file))
                    if "result" not in data_task:
                        raise ValueError("result information is missing in {0}".format(self.task_file))
                    if not os.path.isdir(data_task["result"]):
                        raise ValueError(
                                "result information is not a directory in {0}".format(self.task_file))
                    if "mail" not in data_task:
                        raise ValueError(
                                "mail information is missing in {0}".format(self.task_file))
                    assert(len(data_task) > 0)
                    data_task_ok = data_task
            except IOError as err:
                err.extra_info("Error cannot open {0}".format(self.task_file)) 
            except AssertionError as err:
                err.extra_info("Json {0} is empty".format(self.task_file))
            except TypeError as err:
                err.extra_info("Check information format in {0}".format(self.task_file))
        except:
            self.logger.error("Failed to read {0}".format(self.task_file))
            self.logger.error(sys.exc_info()[1])
            shutil.move(self.task_file, self.error_dir +  
                        os.path.basename(self.task_file))
        return data_task_ok


    def send_fastq(self, history_id, path):
        """Send fastq file
        """
        collection_description = {'collection_type': 'list',
                                   'element_identifiers': [],
                                   'name': "collection_{0}".format(str(os.getpid()))}
        for i,file in enumerate(glob.glob('{0}/*.f*q*'.format(path))):
            dataset = self.gi.tools.upload_file(file, history_id)
            # Add dataset in the collection
            collection_description['element_identifiers'].append(
                {'id': dataset['outputs'][0]["id"],
                'name': "element {0}".format(i),
                'src': 'hda'})
        return collection_description

    def paired_process(self, history):
        """
        """
        dataset_map = {}
        # Upload data
        fasta_dataset_R1 = self.gi.tools.upload_file(
                        self.data_task["contaminant_R1"],
                        history['id'])
        fasta_dataset_R2 = self.gi.tools.upload_file(
            self.data_task["contaminant_R2"],
            history['id'])
        # Upload fastq
        collection_description_R1 = self.send_fastq(
            history['id'], self.data_task["path_R1"])
        collection_description_R2 = self.send_fastq(
            history['id'], self.data_task["path_R2"])
        # Create collection
        collection_R1 = self.gi.histories.create_dataset_collection(
            history['id'], collection_description_R1)
        collection_R2 = self.gi.histories.create_dataset_collection(
            history['id'], collection_description_R2)
        # Get the workflow
        workflow = self.gi.workflows.get_workflows(
            name="masque_paired_end_" + self.data_task["type"])
        detailworkflow = self.gi.workflows.show_workflow(
            workflow[0]['id'])
        return workflow, dataset_map

    def single_process(self, history):
        """
        """
        dataset_map = {}
        # Upload data
        fasta_dataset = self.gi.tools.upload_file(
            self.data_task["contaminant"],
            history['id'])
        # Upload fastq
        collection_description = self.send_fastq(
            history['id'], self.data_task["path"])
        # Create collection
        collection = self.gi.histories.create_dataset_collection(
            history['id'], collection_description)
        # Get the workflow
        workflow = self.gi.workflows.get_workflows(
            name="masque_single_end_" + self.data_task["type"])
        detailworkflow = self.gi.workflows.show_workflow(
            workflow[0]['id'])
        # Get fastq input
        input_collection = self.gi.workflows.get_workflow_inputs(
            workflow[0]['id'], label='Input Dataset Collection')[0]
        # Dataset input
        dataset_map[input_collection] = {'id' : collection['id'],
                                         'src' : 'hdca'}
        # Get contaminant input                                                 
        input_fasta = self.gi.workflows.get_workflow_inputs(
            workflow[0]['id'], label='Input Dataset')[0]
        # Dataset input
        dataset_map[input_fasta] = {'id' : fasta_dataset['outputs'][0]['id'],
                                    'src' : 'hda'}
        return workflow, dataset_map

    def reconnect(self):
        """Reconnect to galaxy
        """
        self.gi = GalaxyInstance(url=self.galaxy_url, key=self.galaxy_key)

    def check_progress(self, history):
        """Check progression
        """
        job_done = False
        try:
            if history:
                # Check status
                while not job_done:
                    print(self.gi.histories.get_status(history['id']))
                    # Success
                    if self.gi.histories.get_status(history['id'])['state'] == "ok":
                        job_done = True
                    # fail 
                    elif self.gi.histories.get_status(history['id'])['state'] == "error": 
                        print(self.gi.histories.get_status(history['id']))
                        break
                    time.sleep(30)                   
        except bioblend.ConnectionError:
            time.sleep(5)
            self.reconnect()
            job_done = self.check_progress(history)
        return job_done

    def send_mail(self, result_file):
        """
        """
        fromaddr = "hub.pasteur@gmail.com"
        toaddr = self.data_task["mail"]
        msg = MIMEMultipart()
        msg['From'] = fromaddr
        msg['To'] = toaddr
        msg['Subject'] = "Shaman result"
        body = "Shaman result is available here"
        msg.attach(MIMEText(body, 'plain')) 
        filename = result_file
        attachment = open(result_file, "rb")
        part = MIMEBase('application', 'octet-stream')
        part.set_payload((attachment).read())
        encoders.encode_base64(part)
        part.add_header('Content-Disposition', "attachment; filename= {0}".format(filename))
        msg.attach(part)
        server = smtplib.SMTP('smtp.gmail.com', 587)
        server.starttls()
        server.login(fromaddr, "c3biRocks")
        text = msg.as_string()
        server.sendmail(fromaddr, toaddr, text)
        server.quit()
    
    def download_result(self, history_id, jeha_id, result_file):
        """
        """
        success = False
        try:
            with open(result_file, 'wb') as result:
                self.gi.histories.download_history(
                    history_id, jeha_id, result)
            if os.stat(result_file).st_size != 0:
                success = True
        except IOError:
            self.logger.error("Error cannot open {0}".format(result_file))
        except requests.exceptions.HTTPError:
            self.logger.error("Error requests try again")
            time.sleep(5)
            success = self.download_result(history_id, jeha_id, result_file)
        return success

    def run(self):
        """Upload, run galaxy workflow and dowload results
        """
        jeha_id = ""
        result_history = None
        dataset_map = None
        data_history_name = 'data_shaman_' + str(os.getpid())+ "_" + str(self.num_job)
        result_history_name = 'shaman_' + str(os.getpid())+ "_" + str(self.num_job)
        # Load json data
        self.data_task = self.load_json()
        if self.data_task:
            # Output result
            result_file = self.data_task["result"] + os.sep + "result.tar.gz"
            # Send data
            try:
                # Create an history
                data_history = self.gi.histories.create_history(
                    name=data_history_name)
                self.logger.info("Load data for {0} : {1}".format(
                    data_history_name, data_history['id']))
                # Send data to the history
                if self.data_task["paired"]:
                    pass
                    #workflow, dataset_map = self.paired_process(data_history)
                    #my_workflow = self.gi.workflows.get_workflows(name="")
                else:
                    workflow, dataset_map = self.single_process(data_history)
                # Parameter
                if self.data_task['host']:
                param = self.gi.workflows.show_workflow(workflow[0]['id'])
                print(param["steps"]["1"])
                # Select genome
                #param["steps"]["1"]["tool_inputs"]["reference_genome"]='{"source": "indexed", "__current_case__": 0, "index":' + self.data_task["host"]+ "}"
                # Start workflow
                #params={"1":param["steps"]["1"]},
            except:
                self.logger.error("Job failed to submit data for history {0}"
                    .format(data_history_name))
                self.logger.error(sys.exc_info()[1])
            # Run workflow
            try:
                if dataset_map:
                    result_history = self.gi.histories.create_history(
                        name=result_history_name)
                    self.logger.info("Load workflow for {0} : {1}".format(
                    data_history_name, result_history['id']))
                    self.gi.workflows.invoke_workflow(
                        workflow[0]['id'],
                        inputs=dataset_map,
                        history_id=result_history['id'])
            except:
                self.logger.error("Job failed at execution for history {0}"
                    .format(result_history_name))
                self.logger.error(sys.exc_info()[1])
                shutil.move(self.task_file, self.error_dir +  
                            os.path.basename(self.task_file))
                # handle error message
                print("Submission failed", file=sys.stderr)
                #delete_history
                #self.gi.histories.delete_history(result_history['id'], purge=True)
            
            if self.check_progress(result_history):
                # Remove reads after success
                self.logger.info("Workflow finished work for {0} : {1}".format(
                    data_history_name, result_history['id']))
                # Build archive
                while jeha_id == "":
                    jeha_id = self.gi.histories.export_history(result_history['id'])
                    time.sleep(5)
                # Download archive
                download_success = self.download_result(result_history['id'], 
                                                       jeha_id, result_file)
                # Send mail
                if os.path.isfile(result_file) and download_success:
                    self.logger.info("Download succeded for {0} : {1}".format(
                    data_history_name, result_history['id']))
                    # Send email
                    #self.send_mail(result_file)
                    # Delete_history
                    #self.gi.histories.delete_history(result_history['id'], purge=True)
                    shutil.move(self.task_file, self.done_dir + 
                            os.path.basename(self.task_file))
                else:
                    self.logger.error("Failed to download result file for {0}"
                        .format(result_history_name))
                    shutil.move(self.task_file, self.error_dir +  
                                os.path.basename(self.task_file))
                    # handle error message
                    print("Job failed during download", file=sys.stderr)
            else:
                self.logger.error("Job is in error state for history {0}"
                    .format(result_history_name))
                shutil.move(self.task_file, self.error_dir +  
                            os.path.basename(self.task_file))
                # handle error message
                print("Job failed during execution", file=sys.stderr)
                #delete_history
                #self.gi.histories.delete_history(result_history['id'], purge=True)

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


def getArguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.add_argument('-u', dest='galaxy_url', type=str, #required=True,
                        default='https://galaxy-dev.web.pasteur.fr',
                        #default='http://127.0.0.1:8080',
                        help='Url to galaxy.')
    parser.add_argument('-k', dest='galaxy_key', type=str, #required=True,
                        default='7ac30484f696937116f960531a05c2b6',
                        #default='a02fb1213a8e22fbfdb4d56e27e41189',
                        help='User galaxy key.')
    parser.add_argument('-w', dest='work_dir', type=isdir, required=True,
                        help='Path to the top directory.')
    args = parser.parse_args()
    return args


def get_log():
    """
    """
    path_log = (tempfile.gettempdir() + os.sep + "shaman_" 
                + str(os.getpid()) + ".log")
    print(path_log)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
    # Create log file
    file_handler = RotatingFileHandler(path_log, 'a', 1000000, 1)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    # Stream in the the console
    ## TO REMOVE IF daemon
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    logger.addHandler(stream_handler)
    return logger


def check_work(todo_dir):
    """Check if a new job need to be done
    """
    return glob.glob('{0}/*.json'.format(todo_dir))

def create_dir(list_dir):
    """
    """
    try:
        for dir_path in list_dir:
            os.mkdir(dir_path)
    except FileExistsError:
        pass


def pandaemonium(logger, galaxy_url, galaxy_key, work_dir):
    """Daemon function that should do something
    """
    todo_dir = work_dir + os.sep + "todo" + os.sep
    doing_dir = work_dir + os.sep + "doing" + os.sep
    done_dir = work_dir + os.sep + "done" + os.sep
    error_dir = work_dir + os.sep + "error" + os.sep
    todo_list = []
    num_job = 0
    # Create important dir
    create_dir([todo_dir, doing_dir, done_dir, error_dir])
    # Start daemon activity
    while True:
        todo_list = check_work(todo_dir)
        if len(todo_list) > 0:
            logger.info("I have a new job todo")
            for task in todo_list: 
                print(task)
                todo_file = doing_dir + os.path.basename(task)
                shutil.move(task, todo_file)  
                djinn = galaxy(logger, todo_file, done_dir, error_dir,
                               galaxy_url, galaxy_key, num_job)
                djinn.start()
                logger.info("task on {0} started".format(todo_file))
                num_job += 1
        time.sleep(5)        
        


def main():
    """Main program
    """
    args = getArguments()
    logger = get_log()
    logger.info("Let's start to work")
    
    #with daemon.DaemonContext():   
    pandaemonium(logger, args.galaxy_url, args.galaxy_key, args.work_dir)


if __name__ == '__main__':
    main()