import os,sys
import argparse
import json
import subprocess
import shutil

parser = argparse.ArgumentParser(description='prepare crab code')
parser.add_argument('-f', dest='file', default='', help='json file input')
parser.add_argument('-m', dest='mode', default='', help='work mode')
args = parser.parse_args()

def get_abbre(name,sample_type,year):
    if sample_type == 'MC':
        return name.split('/')[1] + '_' + year
    elif sample_type == 'data':
        return name.split('/')[1] + '_' + name.split('/')[2].split('-')[0]

def prepare_crab(name,sample_type,year,era,outLFNDirBase):

    abbre_name = get_abbre(name,sample_type,year) 
    if not os.path.exists('crabcode_' + year):
        os.mkdir("crabcode_" + year)

    print ("------> preparing submit code for",abbre_name)
    with open('crabcode_' + year + '/' + abbre_name + '_cfg.py', 'w+') as f:
        f.write('from WMCore.Configuration import Configuration \n\n')

        f.write('config = Configuration()\n')
        f.write('config.section_("General")\n')
        f.write('config.General.requestName = "' + abbre_name + '"\n')
        f.write('config.General.transferLogs = False \n')
        f.write('config.General.workArea = "crab' + year + '"\n\n')

        f.write('config.section_("JobType")\n')
        f.write('config.JobType.pluginName = "Analysis"\n')
        f.write('config.JobType.psetName = "PSet.py"\n')
        f.write('config.JobType.scriptExe = "./go_crab.sh" \n')
        f.write('config.JobType.inputFiles = ["./haddnano.py","./TTC_postproc.py","./keep_and_drop.txt"] #hadd nano will not be needed once nano tools are in cmssw \n')
        f.write('config.JobType.scriptArgs = ["isdata=' + sample_type + '","year=' + year + '","era=' + era + '"] \n')
        f.write('config.JobType.sendPythonFolder  = True\n')
        f.write('config.JobType.allowUndistributedCMSSW = True \n')
        # f.write('config.JobType.maxJobRuntimeMin = 4320 \n\n')

        f.write('config.section_("Data")\n')
        f.write('config.Data.inputDataset = "' + name + '" \n')
        f.write('#config.Data.inputDBS = "phys03"\n')
        f.write('config.Data.inputDBS = "global"\n')
        if sample_type == 'MC':
            f.write('config.Data.splitting = "FileBased"\n')
            f.write('config.Data.unitsPerJob = 1\n')
        elif sample_type == 'data':
            f.write('config.Data.splitting = "LumiBased"\n')
            f.write('config.Data.unitsPerJob = 80\n')
        f.write('#config.Data.splitting = "EventAwareLumiBased" \n')
        # f.write('config.Data.splitting = "Automatic" \n')

        if sample_type == 'MC':
            pass
        elif year == '2018':
            if 'MuonEG_Run2018D' not in abbre_name:
                f.write('config.Data.lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt" \n\n')
            else:
                f.write('config.Data.lumiMask = "" \n\n')
        elif year == '2017':
            f.write('config.Data.lumiMask = "Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt" \n\n')


        f.write('config.Data.outLFNDirBase ="{outLFNDirBase}' + sample_type + '/' + year + '"\n')
        f.write('config.Data.publication = False\n')
        f.write('config.Data.ignoreLocality = True\n')
        f.write('config.Data.allowNonValidInputDataset = True\n')
        f.write('config.Data.outputDatasetTag = "' + abbre_name + '" \n\n')

        f.write('config.section_("Site")\n')
        f.write('config.Site.storageSite = "T2_CN_Beijing"\n')
        f.write('config.Site.whitelist = ["T2_US_MIT","T2_US_Wisconsin","T2_US_Purdue","T2_US_UCSD","T2_US_Caltech"] \n')
        f.close()

def submit(name,sample_type,year):

    abbre_name = get_abbre(name,sample_type,year)

    if not os.path.exists('crabcode_{year}/{abbre_name}_cfg.py'):
        print ("crabcode for ",abbre_name," not existed, \033[31mskipping\033[0m")
        return True

    r=subprocess.run(args="crab submit -c crabcode_{year}/{abbre_name}_cfg.py",shell=True,stdout=subprocess.PIPE,encoding='utf-8')
    if 'Success' in r.stdout:
        print ("--------> submit info:","submit crab jobs for",abbre_name)
    else:
        print ("--------> submit info:","\033[31mfail\033[0m to submit for",abbre_name)

def kill(name,sample_type,year):

    abbre_name = get_abbre(name,sample_type,year)

    if not os.path.exists('crab{year}/crab_{abbre_name}'):
        print ("crab log for ",abbre_name," not existed, \033[31mskipping\033[0m \n")
        return True

    r=subprocess.run(args="crab kill -d crab{year}/crab_{abbre_name}" ,shell=True,stdout=subprocess.PIPE,encoding='utf-8')
    print (r.stdout,'\n')

    shutil.rmtree('crab{year}/crab_{abbre_name}')
    print ('crab{year}/crab_{abbre_name} has been removed')


def status(name,sample_type,year):

    abbre_name = get_abbre(name,sample_type,year)

    if not os.path.exists('crab{year}/crab_{abbre_name}'):
        print ("crab log for ",abbre_name," not existed, \033[31mskipping\033[0m \n")
        return True

    r=subprocess.run(args="crab status -d crab{year}/crab_{abbre_name}" ,shell=True,stdout=subprocess.PIPE,encoding='utf-8')
    print (r.stdout,'\n')


def hadd_help(name,sample_type,year,store_path):

    abbre_name = get_abbre(name,sample_type,year)
    first_name = name.split('/')[1]

    if os.path.exists('{abbre_name}.root'):
        print ('{abbre_name} already existed, \033[31mskipping\033[0m')
        return True

    if not (os.path.exists('{store_path}/{sample_type}/{year}/{first_name}/{abbre_name}')):
        print ('results for {abbre_name} not existed in {store_path}/{sample_type}/{year}/{first_name}/{abbre_name}, \033[31mskipping\033[0m\n')
        return True
    
    if not (len(os.listdir('{store_path}/{sample_type}/{year}/{first_name}/{abbre_name}')) == 1 ):
        print('more than 1 result for {abbre_name}, Please check {store_path}/{sample_type}/{year}/{first_name}/{abbre_name}\n')
        return True

    run_number = os.listdir('{store_path}/{sample_type}/{year}/{first_name}/{abbre_name}')[0]
    path = '{store_path}/{sample_type}/{year}/{first_name}/{abbre_name}/{run_number}/0000/'
    print ('hadding root files in {path}')
    r=subprocess.run(args="haddnano.py {abbre_name}.root {path}/*.root ", shell=True,stdout=subprocess.PIPE,encoding='utf-8')
    
    if os.path.exists('{abbre_name}.root'):
        print ('hadd \033[32mcomplete\033[0m, please check {abbre_name}.root\n')
    else:
        print ('hadd \033[31m fail \033[0m!!')

def report_lumi(name,sample_type,year):

    abbre_name = get_abbre(name,sample_type,year)
    if not os.path.exists('crab{year}/crab_{abbre_name}'):
        print ("crab log for ",abbre_name," not existed, \033[31mskipping\033[0m \n")
        return True

    r=subprocess.run(args="crab report -d crab{year}/crab_{abbre_name}" ,shell=True,stdout=subprocess.PIPE,encoding='utf-8')
    print (r.stdout,'\n')

    if not os.path.exists('lumi_{year}'):
        os.mkdir('lumi_{year}')
    
    shutil.copy('crab{year}/crab_{abbre_name}/results/notFinishedLumis.json', 'lumi_{year}/{abbre_name}.json')

def resubmit(name,sample_type,year):

    abbre_name = get_abbre(name,sample_type,year)

    if not os.path.exists('crab{year}/crab_{abbre_name}'):
        print ("crab log for ",abbre_name," not existed, \033[31mskipping\031[0m \n")
        return True

    print ("resubmitting {abbre_name}\n")
    r = subprocess.run(args="crab resubmit -d crab{year}/crab_{abbre_name}" ,shell=True,stdout=subprocess.PIPE,encoding='utf-8')
    print (r.stdout,"\n")


if __name__=='__main__':
    
    store_path = '/eos/user/s/sdeng/TTC/test/'
    outLFNDirBase = '/store/user/sdeng/TTC/test/'

    with open(args.file, "r") as f:
        jsons = json.load(f)
        f.close()

    if args.mode == 'prepare':
        for dataset in jsons:
            prepare_crab(dataset['name'], dataset['type'], str(dataset['year']), dataset['era'], outLFNDirBase)
    
    if args.mode == 'submit':
        for dataset in jsons:
            submit(dataset['name'], dataset['type'], str(dataset['year']))

    if args.mode == 'kill':
        for dataset in jsons:
            kill(dataset['name'], dataset['type'], str(dataset['year']))
    
    if args.mode == 'status':
        for dataset in jsons:
            status(dataset['name'], dataset['type'], str(dataset['year']))

    if args.mode == 'hadd':
        for dataset in jsons:
            hadd_help(dataset['name'], dataset['type'], str(dataset['year']), store_path)

    if args.mode == 'report':
        for dataset in jsons:
            if dataset['type'] == 'data':
                report_lumi(dataset['name'], dataset['type'], str(dataset['year']))
    
    if args.mode == 'resubmit':
        for dataset in jsons:
            resubmit(dataset['name'], dataset['type'], str(dataset['year']))
