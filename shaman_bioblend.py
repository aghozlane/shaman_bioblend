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
#import keyring
#import tarfile
import zipfile
import lockfile
import datetime
import socket

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

class galaxy(Thread):

    def __init__(self, logger, task_file, doing_dir, done_dir, error_dir,
                 galaxy_url, galaxy_key, num_job, https_mode, delete_mode):
        Thread.__init__(self)
        self.logger = logger
        self.galaxy_url = galaxy_url
        self.galaxy_key = galaxy_key
        self.logger.info("Starting galaxy instance for {0} : {1}".format(
                    galaxy_url, galaxy_key))
        self.gi = GalaxyInstance(url=galaxy_url, key=galaxy_key)
        self.logger.info("{0}".format(self.gi))
        self.logger.info("Connection obtained for {0} : {1}".format(
                    galaxy_url, galaxy_key))
        self.gi.verify = https_mode
        self.task_file = task_file
        self.doing_dir = doing_dir
        self.done_dir = done_dir
        self.error_dir = error_dir
        self.num_job = num_job
        self.dataset_ids = []
        self.delete_mode = delete_mode

    def load_json(self):
        """Load and validate Json
        """
        data_task_ok = None
        try:
            try:
                with open(self.task_file, "rt") as task:
                    data_task = json.load(task)
                    if isinstance(data_task, (list, tuple)):
                        data_task = data_task[0]
                    data_task["name"] = os.path.splitext(os.path.basename(self.task_file))[0]
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
                            if "pattern_R1" not in data_task:
                                raise ValueError(
                                "pattern_R1 information is missing in {0}".format(self.task_file))
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
                            "contaminant information is not a file for {0}".format(self.task_file))
                    if "host" not in data_task:
                        raise ValueError("host information is missing in {0}".format(self.task_file))
                    if "type" not in data_task:
                        raise ValueError("type information is missing in {0}".format(self.task_file))
                    if data_task["type"] not in ["16S", "18S", "23S_28S",
                                                 "ITS", "WGS"]:
                        raise ValueError("type information is incorrect in {0}".format(self.task_file))
                    if "mail" not in data_task:
                        raise ValueError(
                                "mail information is missing in {0}".format(self.task_file))
                    assert(len(data_task) > 0)
                    data_task_ok = data_task
            except IOError as err:
                err.extra_info("Error cannot open {0}".format(self.task_file)) 
            except AssertionError as err:
                err.extra_info("Json {0} is empty".format(self.task_file))
            #except TypeError as err:
            #    err.extra_info("Check information format in {0}".format(self.task_file))
        except:
            self.logger.error("Failed to read {0}".format(self.task_file))
            self.logger.error(sys.exc_info()[1])
            shutil.move(self.task_file, self.error_dir +  
                        os.path.basename(self.task_file))
        return data_task_ok

    def dump_json(self):
        """Dump json file with galaxy info
        """
        todo_file = self.doing_dir + os.path.basename(self.task_file)
        os.remove(self.task_file)  
        try:
            with open(todo_file, "wt") as task:
                json.dump(self.data_task, task)
        except IOError:
            self.logger.error("Failed to write {0}".format(self.task_file))
        self.task_file = todo_file


    def check_file_size(self, path):
        """Check if no file above 2Gb
        """
        large_file_size = False
        for file in glob.glob('{0}/*.f*q*'.format(path)):
            if os.path.getsize(file) > 2000000000:
                large_file_size = True
        return large_file_size

    def send_fastq(self, history_id, path, lib=None):
        """Send fastq file
        """
        collection_description = {'collection_type': 'list',
                                   'element_identifiers': [],
                                   'name': "collection_{0}".format(str(os.getpid()))}
        for i,fastq_file in enumerate(sorted(glob.glob('{0}/*.f*q*'.format(path)))):
            retry = 0
            send_is_ok = False
            while not send_is_ok and retry <= 5 :
                try:
                    if lib:
                        lib_dataset = self.gi.libraries.upload_file_from_local_path(
                                    lib['id'], fastq_file)
                        # move the data in the history
                        dataset = self.gi.histories.upload_dataset_from_library(
                                        history['id'], lib_dataset[0]['id'])
                    else:
                        dataset = self.gi.tools.upload_file(fastq_file, history_id)
                    if 'outputs' in dataset:
                        if "id" in dataset['outputs'][0]:
                            send_is_ok = True
                        else:
                            print("retry")
                            time.sleep(5)
                            retry += 1 
                    else:
                        print("retry")
                        time.sleep(5)
                        retry += 1 
                except:
                    print("retry")
                    time.sleep(5)
                    retry += 1 
            # Add dataset in the collection
            collection_description['element_identifiers'].append(
                {'id': dataset['outputs'][0]["id"],
                'name': "element {0}".format(i),
                'src': 'hda'})
            if self.delete_mode:
                os.remove(fastq_file)
        return collection_description#, i

    def paired_process(self, history, lib=None):
        """
        """
        dataset_map = {}
        # Upload data
        fasta_dataset = self.gi.tools.upload_file(
            self.data_task["contaminant"], history['id'])
        # Upload fastq
        # , count_r1
        collection_description_R1 = self.send_fastq(
            history['id'], self.data_task["path_R1"], lib)
        # , count_r2
        collection_description_R2 = self.send_fastq(
            history['id'], self.data_task["path_R2"], lib)
        # Create collection
        collection_R1 = self.gi.histories.create_dataset_collection(
            history['id'], collection_description_R1)
        collection_R2 = self.gi.histories.create_dataset_collection(
            history['id'], collection_description_R2)
        # Get the workflow
        if self.data_task['host'] != "":
            workflow = self.gi.workflows.get_workflows(
                name="masque_paired_end_" + self.data_task["type"])
        else:
            workflow = self.gi.workflows.get_workflows(
                name="masque_paired_end_" + self.data_task["type"]
                + "_short")
        #detailworkflow = self.gi.workflows.show_workflow(
        #    workflow[0]['id'])
        # Get fastq input
        input_collection_R1 = self.gi.workflows.get_workflow_inputs(
            workflow[0]['id'], label='reads_dataset_collection_R1')[0]
        input_collection_R2 = self.gi.workflows.get_workflow_inputs(
            workflow[0]['id'], label='reads_dataset_collection_R2')[0]
        # Dataset input
        dataset_map[input_collection_R1] = {'id':collection_R1['id'],
                                            'src':'hdca'}
        dataset_map[input_collection_R2] = {'id':collection_R2['id'],
                                            'src':'hdca'}
        # Get contaminant input                                                 
        input_fasta = self.gi.workflows.get_workflow_inputs(
            workflow[0]['id'], label='contaminant_dataset')[0]
        # Dataset input
        dataset_map[input_fasta] = {'id':fasta_dataset['outputs'][0]['id'],
                                    'src':'hda'}
        return workflow, dataset_map#, count_r1 + count_r2

    def single_process(self, history, lib=None):
        """Load fastq, create collection and identify workflow
        """
        dataset_map = {}
        # Upload data
        fasta_dataset = self.gi.tools.upload_file(
            self.data_task["contaminant"], history['id'])
        # Upload fastq
        #, count_fastq
        collection_description = self.send_fastq(history['id'],
                                                 self.data_task["path"],
                                                 lib)
        # Create collection
        collection = self.gi.histories.create_dataset_collection(
            history['id'], collection_description)
        # Get the workflow
        if self.data_task['host'] != "":
            workflow = self.gi.workflows.get_workflows(
                name="masque_single_end_" + self.data_task["type"])
        else:
            workflow = self.gi.workflows.get_workflows(
                name="masque_single_end_" + self.data_task["type"]
                + "_short")
        #detailworkflow = self.gi.workflows.show_workflow(
        #    workflow[0]['id'])
        # Get fastq input
        input_collection = self.gi.workflows.get_workflow_inputs(
            workflow[0]['id'], label='reads_dataset_collection')[0]
        # Dataset input
        dataset_map[input_collection] = {'id' : collection['id'],
                                         'src' : 'hdca'}
        # Get contaminant input                                                 
        input_fasta = self.gi.workflows.get_workflow_inputs(
            workflow[0]['id'], label='contaminant_dataset')[0]
        # Dataset input
        dataset_map[input_fasta] = {'id' : fasta_dataset['outputs'][0]['id'],
                                    'src' : 'hda'}
        return workflow, dataset_map#, count_fastq

    def reconnect(self):
        """Reconnect to galaxy
        """
        self.gi = GalaxyInstance(url=self.galaxy_url, key=self.galaxy_key)
        self.gi.verify = False

    def get_unique(self, seq):
        # Not order preserving
        return {}.fromkeys(seq).keys()

    def check_progress(self, history, glob_progress=0.0, prev_progress=0.0):
        """Check progression
        """
        countdown = 0
        progress_file = (self.doing_dir + os.sep + self.data_task["name"]
                         + "_progress.txt")
        error_file = (self.error_dir + os.sep + self.data_task["name"]
                         + "_error.txt")
        job_done = False
        error_list = []
        error_mess_list = []
        try:
            # Check status
            while not job_done:
                progress_story = self.gi.histories.get_status(history['id'])
                #new_progress = float(self.gi.histories.get_status(history['id'])['percent_complete'])
                new_progress = float(progress_story['percent_complete'])
                if prev_progress > new_progress:
                    glob_progress = 200.0 + new_progress
                    prev_progress = glob_progress
                else:
                    glob_progress = glob_progress + (new_progress - prev_progress)
                    prev_progress = new_progress
                with open(progress_file, "wt") as progress:
                    progress.write("{0}".format(glob_progress / 3.0))
                # print("progression {0}".format(glob_progress))
                # print("progression {0}".format(glob_progress / 3.0))
                
                # Success
                if progress_story['state'] == "ok" and (glob_progress == 100.0 or glob_progress == 300.0):
                    job_done = True
                    self.logger.info(progress_story)
                elif progress_story['state'] == "ok" and countdown >= 10:
                    job_done = True
                    self.logger.info(progress_story)
                elif progress_story['state'] == "ok":
                    countdown += 1
                # fail 
                elif progress_story['state'] == "error" or progress_story['state_details']['error'] > 0: 
                    self.logger.error(progress_story)
                    error_datasets = self.gi.histories.show_history(history['id'])['state_ids']['error']
                    error_jobid = [self.gi.histories.show_dataset_provenance(history['id'], dataset_id)['job_id']
                                   for dataset_id in error_datasets]
                    # Unique jobid in error
                    error_jobid.sort()
                    error_jobid = self.get_unique(error_jobid)
                    # Get error message
                    for job_id in error_jobid:
                        error_list.append(self.gi.jobs.show_job(job_id, full_details=True))
                    # Write error message
                    with open(error_file, "wt") as error_log:
                        for error in error_list:
                            error_mess = ("tool_id:{1}{0}error:{0}{2}"
                                            .format(os.linesep, 
                                                    error['tool_id'],
                                                    error['stderr']))
                            error_log.write(error_mess)
                            error_mess_list.append(error_mess)

                    message = ("The workflow failed during progression for the "
                               "key {0}.{1}{2}"
                               .format(self.data_task["name"].replace("file", ""),
                                os.linesep, os.linesep.join(error_mess_list)))
                    # error mail
                    self.send_mail(message)
                    break
                else:
                    self.logger.info(progress_story)
                time.sleep(10)                   
        except bioblend.ConnectionError:
            time.sleep(5)
            self.reconnect()
            job_done = self.check_progress(history, glob_progress, prev_progress)
        except IOError:
            self.logger.error("Error cannot open {0} or {1}"
                              .format(progress_file, error_file))
            job_done = self.check_progress(history, glob_progress, prev_progress)
        return job_done

    # def get_members(self, tar, prefix):
    #     """Identify element in the tar.gz
    #     """
    #     if not prefix.endswith(os.sep):
    #         prefix += os.sep
    #     offset = len(prefix)
    #     for tarinfo in tar.getmembers():
    #         if tarinfo.name.startswith(prefix):
    #             tarinfo.name = tarinfo.name[offset:]
    #             yield tarinfo

    # def zipdir(self, path, ziph):
    #     """List elements in the directory and write each file
    #     """
    #     # ziph is zipfile handle
    #     for root, dirs, files in os.walk(path):
    #         for file in files:
    #             ziph.write(os.path.join(root, file), file)

    def zip_archive(self, list_downloaded_files, zip_file):
        """Extract tar file and build a clean zip file
        """
        #if not os.path.isdir(path):
        #    os.mkdir(path)
        #with tarfile.open(result_file) as tar:
        #    tar.extractall(path, self.get_members(tar, "datasets"))
        #print("remove " + result_file)
        #os.remove(result_file)
        # Rename
        # for file in glob.glob('{0}/*'.format(result_dir)):
        #     split_name = os.path.splitext(os.path.basename(file))
        #     if split_name[1] == ".tabular":
        #         shutil.move(file,  path + os.sep + 
        #                     "_".join(split_name[0].split("_")[0:3]) + ".tsv")
        #     elif split_name[1] == ".biom1":
        #         shutil.move(file,  path + os.sep + 
        #                     "_".join(split_name[0].split("_")[0:2]) + ".biom")
        #     else:
        #         resplit = split_name[0].split("_")
        #         shutil.move(file, path + os.sep +
        #             "_".join(resplit[0:len(resplit)-1]) + split_name[1])
        #try:
        with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for file in list_downloaded_files:
                zipf.write(file, os.path.basename(file))
        #except IOError:
        #    self.logger.error("Error cannot open {0}".format(zip_file))

    def send_mail(self, message, result_file=None):
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
        msg.attach(MIMEText(message, 'plain'))
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
        #smtplib.SMTPSenderRefused:
        #except IOError:
        #    self.logger.error("Error cannot open {0}".format(result_file))
    
    # def download_result(self, history_id, jeha_id, result_file):
    #     """Download tar archive from galaxy
    #     """
    #     success = False
    #     try:
    #         with open(result_file, 'wb') as result:
    #             self.gi.histories.download_history(
    #                 history_id, jeha_id, result)
    #         if os.stat(result_file).st_size != 0:
    #             success = True
    #     except IOError:
    #         self.logger.error("Error cannot open {0}".format(result_file))
    #     except requests.exceptions.HTTPError:
    #         self.logger.error("Error requests try again")
    #         time.sleep(5)
    #         success = self.download_result(history_id, jeha_id, result_file)
    #     return success

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
            self.logger.error("Failed to download {0}".format(res))
            success = False
        except AssertionError:
            self.logger.error("Failed to download {0}, no file matched".format(res))
            success = False
        return success, list_downloaded_files

    def run(self):
        """Upload, run galaxy workflow and dowload results
        """
        list_result = {}
        #count_fastq = 0
        lib = None
        param = None
        jeha_id = ""
        result_history = None
        dataset_map = None
        data_history_name = ('data_shaman_' + str(os.getpid())+ "_" + 
                             str(self.num_job))
        lib_name = 'lib_shaman_' + str(os.getpid())+ "_" + str(self.num_job)
        result_history_name = ('shaman_' + str(os.getpid())+ "_" + 
                                str(self.num_job))
        # Load json data
        self.logger.info("Start reading {0}".format(
                    self.task_file))
        self.data_task = self.load_json()
        self.logger.info("Done reading {0}".format(
                    self.task_file))

        # Add galaxy info
        self.data_task['data_history_name'] = data_history_name
        self.data_task['result_history_name'] = result_history_name
        self.logger.info("Starting dump of {0}".format(
                    self.task_file))
        self.dump_json()
        self.logger.info("Done dumping of {0}".format(
                    self.task_file))
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
        # Check file size
        if self.data_task:
            # Output result
            #result_file = self.done_dir + os.sep + self.data_task["name"] + ".tar.gz"
            
            # Send data
            try:
                # Create an history
                self.logger.info("Starting new history {0}".format(
                    data_history_name))
                data_history = self.gi.histories.create_history(
                    name=data_history_name)
                self.logger.info("Load data for {0} : {1}".format(
                    data_history_name, data_history['id']))
                # Send data to the history
                if self.data_task["paired"]:
                    if (self.check_file_size(self.data_task["path_R1"]) or 
                        self.check_file_size(self.data_task["path_R2"])):
                        lib = self.gi.libraries.create_library(lib_name)
                        #, count_fastq
                        workflow, dataset_map = self.paired_process(
                            data_history, lib)
                    else:
                        #, count_fastq
                        workflow, dataset_map = self.paired_process(data_history)
                else:
                    # Check file size
                    if self.check_file_size(self.data_task["path"]):
                        lib = self.gi.libraries.create_library(lib_name)
                        #, count_fastq
                        workflow, dataset_map = self.single_process(
                            data_history, lib)
                    else:
                        #, count_fastq
                        workflow, dataset_map = self.single_process(data_history)
            except:
                self.logger.error("Shaman failed to submit data for the history {0}"
                    .format(data_history_name))
                self.logger.error(sys.exc_info()[1])
                shutil.move(self.task_file, self.error_dir +  
                            os.path.basename(self.task_file))
                message = "Shaman failed to submit data for the history {0}".format(self.data_task["name"].replace("file", ""))
                self.send_mail(message)
                # if lib:
                #    self.gi.libraries.delete_library(lib['id'])
                # #delete_history
                # if data_history:
                #     self.gi.histories.delete_history(data_history['id'], purge=True)
                data_history = None
            # Run workflow
            if data_history:
                if self.check_progress(data_history):
                    try:
                        if dataset_map:
                            #result_history = data_history
                            result_history = self.gi.histories.create_history(
                                                 name=result_history_name)
                            self.logger.info("Load workflow for {0} : {1}".format(
                            data_history_name, result_history['id']))
                            align_dict = {
                            'id':self.data_task["aKmin"],
                            'strand':self.data_task["annotationstrand"]
                            }
                            clustering_dict = {
                                'id':self.data_task["clusteringthreshold"],
                                'strand':self.data_task["clusteringstrand"]
                            }
                            annot_dict = {
                                'aKmin':self.data_task["aKmin"],
                                'aPmin':self.data_task["aPmin"],
                                'aPmax':self.data_task["aPmax"],
                                'aCmin':self.data_task["aCmin"],
                                'aCmax':self.data_task["aCmax"],
                                'aOmin':self.data_task["aOmin"],
                                'aOmax':self.data_task["aOmax"],
                                'aFmin':self.data_task["aFmin"],
                                'aFmax':self.data_task["aFmax"],
                                'aGmin':self.data_task["aGmin"],
                                'aGmax':self.data_task["aGmax"],
                                'aSmin':self.data_task["aSmin"]
                            }
                            quality_dict = {
                                'q': self.data_task["phredthres"],
                                'p': self.data_task["mincorrect"],
                                'l': self.data_task["minreadlength"]
                            }
                            derep_dict = {
                                'derep_method': self.data_task["dreptype"],
                                'minseqlength': self.data_task["minampliconlength"]
                            }
                            # paired end with host
                            if self.data_task['host'] != "" and self.data_task["paired"]:
                                params = {
                                    "3":{"reference_genome|index":self.data_task["host"]},
                                    "5": quality_dict,
                                    "7":{
                                    'pattern|sub_pattern': self.data_task["pattern_R1"],
                                    'max_amplicon_length':self.data_task["maxampliconlength"]},
                                    "9": derep_dict,
                                    "10":{'sorting_mode|minsize':self.data_task["minabundance"]},
                                    #Clustering
                                    "12":clustering_dict
                                         }
                                if self.data_task["type"] == "16S" :
                                    params.update(
                                        {
                                        #Count matrix
                                        "16":clustering_dict,
                                        #Greengenes
                                        "14":align_dict,
                                        "18":annot_dict,
                                        #Silva
                                        "15":align_dict,
                                        "19":annot_dict,
                                        #Extract Result
                                        "23":{'paired|pattern': self.data_task["pattern_R1"]}
                                        }
                                    )
                                elif self.data_task["type"] == "18S" or self.data_task["type"] == "23S_28S":
                                    params.update(
                                        {
                                        # Count matrix
                                        "16":clustering_dict,
                                        # Silva
                                        "14":align_dict,
                                        "17":annot_dict,
                                        # Extract Result
                                        "20":{'paired|pattern': self.data_task["pattern_R1"]}
                                        }
                                    )
                                else:
                                    params.update(
                                        {
                                        # Count matrix
                                        "17":clustering_dict,
                                        # Findley
                                        "14":align_dict,
                                        "19":annot_dict,
                                        # Underhill
                                        "15":align_dict,
                                        "20":annot_dict,
                                        # Unite
                                        "18":align_dict,
                                        "21":annot_dict,
                                        # Extract Result
                                        "28":{'paired|pattern': self.data_task["pattern_R1"]}
                                        }
                                    )
                                self.gi.workflows.invoke_workflow(
                                    workflow[0]['id'], inputs=dataset_map,
                                    params=params, history_id=result_history['id'])
                            # paired end no host
                            elif self.data_task['host'] == "" and self.data_task["paired"]:
                                params={
                                    "3":quality_dict,
                                    "5":{'pattern|sub_pattern': self.data_task["pattern_R1"]},
                                    "7": derep_dict,
                                    "8":{'sorting_mode|minsize':self.data_task["minabundance"]},
                                    #Clustering
                                    "10":clustering_dict,
                                    }
                                if self.data_task["type"] == "16S":
                                    params.update(
                                        {
                                        #Count matrix
                                        "14":clustering_dict,
                                        #Greengenes
                                        "12":align_dict, "16":annot_dict,
                                        #Silva
                                        "13":align_dict, "17":annot_dict,
                                        #Extract result
                                        "21":{'paired|pattern' : self.data_task["pattern_R1"]}
                                        }
                                    )
                                elif self.data_task["type"] == "18S":
                                    params.update(
                                        {
                                        # Count matrix
                                        "14":clustering_dict,
                                        # Silva
                                        "12":align_dict, "15":annot_dict,
                                        # Extract result
                                        "17":{'paired|pattern': self.data_task["pattern_R1"]}
                                        }
                                    )
                                elif self.data_task["type"] == "23S_28S":
                                    params.update(
                                        {
                                        # Count matrix
                                        "14":clustering_dict,
                                        # Silva
                                        "12":align_dict, "15":annot_dict,
                                        # Extract result
                                        "18":{'paired|pattern': self.data_task["pattern_R1"]}
                                        }
                                    )
                                else:
                                    params.update(
                                        {
                                        # Count matrix
                                        "15":clustering_dict,
                                        # Findley
                                        "12":align_dict, "17":annot_dict,
                                        # Underhill
                                        "13":align_dict, "18":annot_dict,
                                        # Unite
                                        "16":align_dict, "19":annot_dict,
                                        # Extract result
                                        "26":{'paired|pattern': self.data_task["pattern_R1"]}
                                        }
                                    )
                                self.gi.workflows.invoke_workflow(
                                    workflow[0]['id'], inputs=dataset_map,
                                    params=params, history_id=result_history['id'])
                            # single end with host
                            elif self.data_task['host'] != "" and not self.data_task["paired"]:
                                params={
                                    "2":{"reference_genome|index":self.data_task["host"]},
                                    "4":quality_dict,
                                    "5":{'max_amplicon_length':self.data_task["maxampliconlength"]},
                                    "7": derep_dict,
                                    "8":{'sorting_mode|minsize':self.data_task["minabundance"]},
                                    #Clustering
                                    "10":clustering_dict,
                                    }
                                if self.data_task["type"] == "16S":
                                    params.update(
                                        {
                                        #Count matrix
                                        "15":clustering_dict,
                                        #Greengenes
                                        "12":align_dict, "16":annot_dict,
                                        #Silva
                                        "13":align_dict, "17":annot_dict
                                        }
                                    )
                                elif self.data_task["type"] == "18S" or self.data_task["type"] == "23S_28S":
                                    params.update(
                                        {
                                        # Count matrix
                                        "14":clustering_dict,
                                        # Silva
                                        "12":align_dict, "15":annot_dict
                                        }
                                    )
                                else:
                                    params.update(
                                        {
                                        # Count matrix
                                        "15":clustering_dict,
                                        # Findley
                                        "12":align_dict, "17":annot_dict,
                                        # Underhill
                                        "13":align_dict, "18":annot_dict,
                                        # Unite
                                        "16":align_dict, "19":annot_dict
                                        }
                                    )
                                self.gi.workflows.invoke_workflow(
                                    workflow[0]['id'], inputs=dataset_map,
                                    params=params,
                                    history_id=result_history['id'])
                            # single end no host
                            else:
                                params={
                                    "2":quality_dict,
                                    "3":{'max_amplicon_length':self.data_task["maxampliconlength"]},
                                    "5":derep_dict,
                                    "6":{'sorting_mode|minsize':self.data_task["minabundance"]},
                                    # Clustering
                                    "8":clustering_dict
                                }
                                if self.data_task["type"] == "16S":
                                    params.update(
                                        {
                                        # Count matrix
                                        "13":clustering_dict,                                        
                                        # Greengenes
                                        "10":align_dict, "14":annot_dict,
                                        # Silva
                                        "11":align_dict, "15":annot_dict
                                        }
                                    )
                                elif self.data_task["type"] == "18S" or self.data_task["type"] == "23S_28S":
                                    params.update(
                                        {
                                        # Count matrix
                                        "12":clustering_dict,
                                        # Silva
                                        "10":align_dict, "13":annot_dict
                                        }
                                    )
                                else:
                                    params.update(
                                        {
                                        # Count matrix
                                        "13":clustering_dict,
                                        # Findley
                                        "10":align_dict, "15":annot_dict,
                                        # Underhill
                                        "11":align_dict, "16":annot_dict,
                                        # Unite
                                        "14":align_dict, "17":annot_dict
                                        }
                                    )
                                self.gi.workflows.invoke_workflow(
                                    workflow[0]['id'], inputs=dataset_map,
                                    params=params,
                                    history_id=result_history['id'])
                    except:
                        self.logger.error("Job failed at execution for the history: {0}"
                            .format(result_history_name))
                        self.logger.error(sys.exc_info()[1])
                        shutil.move(self.task_file, self.error_dir +  
                                    os.path.basename(self.task_file))
                        message = "Workflow failed to start for the key: {0}".format(self.data_task["name"].replace("file", ""))
                        self.send_mail(message)
                        # # Delete history
                        # if lib:     
                        #    self.gi.libraries.delete_library(lib['id'])
                        # #delete_history
                        # if data_history:
                        #     self.gi.histories.delete_history(data_history['id'], purge=True)
                        # if result_history:
                        #     self.gi.histories.delete_history(result_history['id'], purge=True)
                        # result_history = None
            #, count_fastq
            if result_history:
                if self.check_progress(result_history, 100.0):
                    # Remove reads after success
                    self.logger.info("Workflow finished work for {0} : {1}".format(
                        data_history_name, result_history['id']))
                    # Download results
                    download_success, list_downloaded_files = self.download_result(
                                            result_history['id'], list_result,
                                            result_dir)
                    # Build archive
                    #jeha_id = self.gi.histories.export_history(result_history['id'], wait=True)

                    # Download archive
                    #download_success = self.download_result(
                    #                        result_history['id'], 
                    #                        jeha_id, result_file)
                    # Send mail
                    #if os.path.isfile(result_file) and download_success:
                    if download_success:
                        self.logger.info("Download succeded for {0} : {1}".format(
                        data_history_name, result_history['id']))
                        # Prepare zip
                        self.zip_archive(list_downloaded_files, zip_file)
                        # Send email
                        # solve file size problem
                        message = ("Shaman result is available for the key {0}"
                            .format(self.data_task["name"].replace("file", "")))
                        self.send_mail(message, zip_file)
                        shutil.move(self.task_file, self.done_dir + 
                                    os.path.basename(self.task_file))
                        # Delete_history
                        if lib:
                           self.gi.libraries.delete_library(lib['id'])
                        self.gi.histories.delete_history(data_history['id'], purge=True)
                        self.gi.histories.delete_history(result_history['id'], purge=True)
                    else:
                        self.logger.error("Failed to download result file for {0}"
                            .format(result_history_name))
                        shutil.move(self.task_file, self.error_dir +  
                                    os.path.basename(self.task_file))
                        message = ("Workflow failed to download the results for the key {0}"
                                   .format(self.data_task["name"].replace("file", "")))
                        self.send_mail(message)
                        # handle error message
                        #print("Job failed during download", file=sys.stderr)
                else:
                    self.logger.error("Workflow failed during progression for the key {0}"
                        .format(result_history_name))
                    #message = "Workflow failed during progression for the key {0}".format(self.data_task["name"].replace("file", ""))
                    #self.send_mail(message)
                    if(os.path.isfile(self.task_file)):
                        shutil.move(self.task_file, self.error_dir +  
                                    os.path.basename(self.task_file))
                    # # delete_library
                    # if lib:
                    #   self.gi.libraries.delete_library(lib['id'])
                    # # delete_history
                    # if data_history:
                    #    self.gi.histories.delete_history(data_history['id'], purge=True)
                    # if result_history:
                    #    self.gi.histories.delete_history(result_history['id'], purge=True)

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
                        default='https://galaxy.pasteur.fr',
                        #default='https://galaxy-dev.web.pasteur.fr',
                        #default='http://127.0.0.1:8080',
                        help='Url to galaxy.')
    parser.add_argument('-k', dest='galaxy_key', type=str, #required=True,
                        default='31f05d9edaa2228b66c538f43b0d5d52',
                        #default=keyring.get_password("galaxy", "aghozlan"),
                        #default='7ac30484f696937116f960531a05c2b6',
                        #default='f293dce7785a77c338db9e8b8df9922c',
                        help='User galaxy key.')
    parser.add_argument('-w', dest='work_dir', type=isdir, required=True,
                        action=FullPaths, help='Path to the top directory.')
    parser.add_argument('-i', dest='interactive_mode', action='store_false',
                        default=True, help='Do not detatch process.')
    parser.add_argument('-s', dest='https_mode', action='store_true',
                        default=False, help='Activate https verification.')
    parser.add_argument('-d', dest='delete_mode', action='store_true',
                        default=False, help='Delete reads provided as input.')
    args = parser.parse_args()
    return args


def get_log(path_log):
    """
    """
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


def pandaemonium(path_log, galaxy_url, galaxy_key, work_dir, https_mode, 
                 delete_mode):
    """Daemon function that should do something
    """
    todo_dir = work_dir + os.sep + "todo" + os.sep
    doing_dir = work_dir + os.sep + "doing" + os.sep
    done_dir = work_dir + os.sep + "done" + os.sep
    error_dir = work_dir + os.sep + "error" + os.sep
    todo_list = []
    num_job = 0

    logger = get_log(path_log)
    logger.info("Let's start to work")
    # Create important dir
    create_dir([todo_dir, doing_dir, done_dir, error_dir])
    # Start daemon activity
    while True:
        todo_list = check_work(todo_dir)
        if len(todo_list) > 0:
            logger.info("I have a new job todo")
            for task in todo_list: 
                #todo_file = doing_dir + os.path.basename(task)
                #shutil.move(task, todo_file)  
                djinn = galaxy(logger, task, doing_dir, done_dir,
                               error_dir, galaxy_url, galaxy_key, num_job,
                               https_mode, delete_mode)
                djinn.start()
                logger.info("task on {0} started".format(task))
                num_job += 1
        time.sleep(5)        
        


def main():
    """Main program
    """
    args = getArguments()
    now = datetime.datetime.now()
    # Check if root
    if os.geteuid() != 0:
        path_log = (tempfile.gettempdir() + os.sep + "shaman_bioblend_" +
                    now.strftime("%Y%m%d_%H%M") + ".log")
    else:
        create_dir(["/var/log/shaman_bioblend"])
        path_log = ("/var/log/shaman_bioblend/shaman_bioblend_" +
                    now.strftime("%Y%m%d_%H%M") + ".log")
    # Start Daemon
    with daemon.DaemonContext(working_directory=os.curdir,
                              pidfile=lockfile.FileLock(path_log),
                              stdout=sys.stdout, stderr=sys.stderr,
                              detach_process=args.interactive_mode):
        print("PID: {0}".format(os.getpid()))
        print("Path to log file: {0}".format(path_log))
        pandaemonium(path_log, args.galaxy_url, args.galaxy_key, args.work_dir,
                     args.https_mode, args.delete_mode)


if __name__ == '__main__':
    main()
