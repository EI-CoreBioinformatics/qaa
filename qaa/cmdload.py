def loadPreCmd(command, is_dependency=True):
    '''
    Used to prefix a shell command that utilises some external software with another command used to load that software
    '''
    if command:
        cc = command.strip()
        if cc != "":
            if is_dependency:
                return "set +u && {} &&".format(cc)
            else:
                return " {} ".format(cc)
    return ""
