def getRemoteData():
    """
    Gets the remote data using remoteData AIAManager class
    """
    from Solar.remoteData import SDOAIAManager
    from astropy import units as u

    sdoaia = SDOAIAManager(aia=SDOAIAManager(times=("2020/5/27",
                                                    "2020/5/28 14:00"),
                                             cadence=1 * u.minute,
                                             aiaPath=""))
    sdoaia.downloadData(force=True)
