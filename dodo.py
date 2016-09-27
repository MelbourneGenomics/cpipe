from tasks.manual_assets import *
from tasks.nectar_assets import *

import os

def task_install():
    swift_credentials = {
        'OS_AUTH_URL',
        'OS_TENANT_ID'
        'OS_TENANT_NAME',
        'OS_PROJECT_NAME',
        'OS_USERNAME',
        'OS_PASSWORD'
    }

    # If the user has all the necessary environment variables set, let them download the cached assets
    # from the object store
    if swift_credentials.issubset(os.environ.keys()):
        task = 'nectar_assets'
    else:
        task = 'manual_assets'

    return {
        'actions': [None],
        'taskdeps': [task]
    }

