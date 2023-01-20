def get_mates(sample):
    """Get number of mates"""
    try:
        if config[sample]['mate2']:
            return "pe"
        else:
            return "se"
    except KeyError:
        return "se"

def get_mates_number(sample):
    """Get library type"""
    try:
        if config[sample]['mate2']:
            return '2'
        else:
            return '1'
    except KeyError:
        return '1'

def get_library_type(sample, sense):
    try:
        if config[sample]['mate2']:
            if str(int(sense)) == '2':
                return "ISR"
            if str(int(sense)) == '1':
                return "ISF"
        else:
            if str(int(sense)) == '2':
                return "SR"
            if str(int(sense)) == '1':
                return "SF"
    except KeyError:
        if str(int(sense)) == '2':
            return "SR"
        if str(int(sense)) == '1':
            return "SF"

def get_mismatches(sample):
    """Get number of mates"""
    try:
        if config[sample]['mismatches']:
            return config[sample]['mismatches']
        else:
            return 0.1
    except KeyError:
        return 0.1
