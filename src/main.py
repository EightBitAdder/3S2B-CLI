from utils import compare, compareAll
from db_utils import download, searchAndFetch, viewIdxTable
from menu.menu_utils import CompareHandler, CompareAllHandler, DownloadHandler, SearchAndFetchHandler, ViewIdxTableHandler
from menu.menu import Menu


def main():
    
    menu = Menu("CRAFTS Lab MFMC Menu")
    
    menu.addMenuItem("V", "View Index Table", ViewIdxTableHandler(viewIdxTable))
    menu.addMenuItem("F", "Search and Fetch", SearchAndFetchHandler(searchAndFetch))
    menu.addMenuItem("D", "Download", DownloadHandler(download))
    menu.addMenuItem("C", "Compare", CompareHandler(compare))
    menu.addMenuItem("A", "Compare All", CompareAllHandler(compareAll))

    menu.run()


if __name__ == "__main__":
    
    main()