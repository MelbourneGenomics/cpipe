def task_nectar_install_dependencies():
    return {
        'actions': ['apt-get install -y gfortran'],
        'uptodate': [False]
    }
