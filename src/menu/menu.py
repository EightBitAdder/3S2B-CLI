import os
import sys
import inspect
from typing import Any, Callable
from prettytable import PrettyTable


class Handler:

    def __init__(self, func: Callable[[Any], Any]) -> None:

        self.func = func


    def handle(self, *args) -> None:

        result = self.func(*args)

        if (result is not None):
            
            print(result)
        


class MenuItem():

    def __init__(self, alias: str, handler: Handler) -> None:

        self.alias = alias
        self.handler = handler
    

class Menu():

    def __init__(self, title: str="") -> None:

        self._title = title
        self._menuItems = {}
        self._addQuit()


    def addMenuItem(self, key: str, alias: str, handler: Handler) -> None:

        self._menuItems[key] = MenuItem(alias, handler)


    def _addQuit(self) -> None:

        def quit() -> None:

            self._clear()
            sys.exit()

        self.addMenuItem("Q", "Quit", Handler(quit))

    
    def _choose(self, key: str) -> None:
        
        if key in self._menuItems:

            menuItem = self._menuItems[key]
            signature = inspect.signature(menuItem.handler.func)

            if (bool(signature.parameters)):

                # TODO:
                # Obviously this is awfuly precarious...
                args = input(f"Enter arg for <<< {menuItem.alias} >>>: ").split(" && ")

                self._clear()
                
                menuItem.handler.handle(*args)
            else:

                menuItem.handler.handle()
        else:

            print("<!> Invalid choice. Please try again.")


    def run(self) -> None:

        while True:

            self._clear()
            self._display()

            choice = input(">>> ").strip().upper()
            
            self._clear()

            try:

                self._choose(choice)
            except Exception as e:

                print(f"<!> Uh oh! There was an error:\n{e}")

            input("\nPress any key to continue...")


    def _clear(self) -> None:
        
        os.system("cls" if os.name == "nt" else "clear")


    def _display(self) -> None:
        
        if (self._title != ""):

            titleTable = PrettyTable()
            titleTable.header = False
            titleTable.junction_char = "-"
            titleTable.add_row([self._title])

            print(titleTable, "\n")
        
        menuTable = PrettyTable()
        menuTable.header = False
        menuTable.junction_char = "-"

        for key, menuItem in self._menuItems.items():

            menuTable.add_row([key, menuItem.alias])

        print(menuTable)