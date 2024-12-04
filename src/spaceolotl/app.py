from shiny import App
from spaceolotl.mod.ui import app_ui as _ui
from spaceolotl.mod.server import server as _server

app = App(_ui, _server)