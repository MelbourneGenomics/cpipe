import platform
import subprocess

def task_ubuntu_dependencies():
    return {
        'actions': [
            '''
               #Fix for docker apt-get
                apt-get update
                apt-get install -y software-properties-common

                # Install gradle repos
                add-apt-repository -y ppa:cwchien/gradle
                apt-get update

                # Install apt-getable things
                apt-get install -y git make poppler-utils zlib1g-dev ncurses-dev gcc g++ patch libssl-dev unzip maven\
                 gradle libcurl4-openssl-dev texinfo openjdk-8-jdk python python-pip mysql-client xorg-dev libreadline-dev libbz2-dev liblzma-dev\
                 libpcre3-dev libsqlite3-dev cpanminus wget curl

            '''],
        'uptodate': [False]
    }


def task_manual_install_dependencies():
    
    if len(subprocess.check_output('cat /etc/*release | egrep "Ubuntu 16"', shell=True)) > 1:
        task = 'ubuntu_dependencies'
    else:
        raise Exception(
            'Unsupported version for manual install. Install the dependencies yourself or obtain NECTAR credentials')

    return {
        'actions': [None],
        'task_dep': [task]
    }
