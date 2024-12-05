from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("amex_merfish_spaceolotl")
except PackageNotFoundError:
    __version__ = "0.0.0_fallback"