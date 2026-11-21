import logging
import tomllib
from pathlib import Path

from config import LoggingConfig

LoggingConfig().initialize()

logger = logging.getLogger(__name__)


class EngineGeometry:
    def __init__(self):
        self.geometry = None

    def from_file(self, path: Path):
        try:
            with open(path, "rb") as f:
                self.geometry = tomllib.load(f)
        except FileNotFoundError:
            logger.warning(f"Geometry file not found: {path}")

        return self

    def get(self):
        return self.geometry


class EngineConfig:
    def __init__(self):
        self.config = None

    def from_file(self, path: Path):
        try:
            with open(path, "rb") as f:
                self.config = tomllib.load(f)
        except FileNotFoundError:
            logger.warning(f"Config file not found: {path}")

        return self

    def get(self):
        return self.config


class EngineSimulation:
    def __init__(self):
        self.config = None
        self.geom = None

    def load_config_from_file(self, config_path: Path, geom_path: Path):
        self.config = EngineConfig().from_file(config_path).get()
        self.geom = EngineGeometry().from_file(geom_path).get()

        return self

    def get_inputs(self):
        if self.config is None:
            logger.warning("No config loaded.")
        else:
            logger.info(f"Loaded config: {self.config}")

        if self.geom is None:
            logger.warning("No geometry loaded.")
        else:
            logger.info(f"Loaded geometry: {self.geom}")


def main():
    config_path = Path("./data/config.toml")
    geom_path = Path("./data/geom.toml")

    sim = EngineSimulation().load_config_from_file(config_path, geom_path)
    sim.get_inputs()


if __name__ == "__main__":
    main()
