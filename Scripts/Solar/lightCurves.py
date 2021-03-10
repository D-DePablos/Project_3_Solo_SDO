import os
import sys

module_path = os.path.abspath(os.path.join(f"{os.getcwd()}/../"))

for module in [module_path]:
    if module not in sys.path:
        print(type(sys.path))
        print(module)
        sys.path.append(module)


def getRemoteData():
    """
    Gets the remote data using remoteData AIAManager class
    """
    from remoteData import SDOAIAManager
    from astropy import units as u

    sdoaia = SDOAIAManager(times=("2020/5/27", "2020/5/28 14:00"),
                           cadence=1 * u.minute,
                           aiaPath="/disk/solar18/ddp/PhD/3_SolO_SDO/unsafe/remoteData/AIA/")
    sdoaia.downloadData(force=True)


getRemoteData()
