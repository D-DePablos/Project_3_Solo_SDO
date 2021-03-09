def getRemoteData():
    """
    Gets the remote data using remoteData AIAManager class
    """
    from .remoteData import SDOAIAManager
    from astropy import units as u

    sdoaia = SDOAIAManager(aia=SDOAIAManager(times=("2020/5/27",
                                                    "2020/5/28 14:00"),
                                             cadence=1 * u.minute,
                                             aiaPath="/disk/solar18/ddp/PhD/Project_3_Solo_SDO/unsafe/remoteData/AIA"))
    sdoaia.downloadData(force=True)


if __name__ == "__main__":
    getRemoteData()
