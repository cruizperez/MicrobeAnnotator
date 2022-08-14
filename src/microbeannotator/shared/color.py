from enum import Enum
from typing import List


class MicrobeAnnotatorColorPalette(Enum):
    """
    Contains the color palette used throughout MicrobeAnnotator.
    """
    MAIN_BLUE = "#003f5c"
    MAIN_PURPLE = "#58508d"
    MAIN_PINK = "#bc5090"
    MAIN_SALMON = "#ff6361"
    MAIN_YELLOW = "#ffa600"
    SECONDARY_ORANGE = "#c7522a"
    SECONDARY_BEIGE = "#e5c185"
    SECONDARY_YELLOW = "#fbf2c4"
    SECONDARY_GREEN = "#74a892"
    SECONDARY_AQUAMARINE = "#008585"
    COMPLEMENTARY_1 = "#f0ead2"
    COMPLEMENTARY_2 = "#f0ead2"
    COMPLEMENTARY_3 = "#f0ead2"
    COMPLEMENTARY_4 = "#dde5b4"
    COMPLEMENTARY_5 = "#adc178"
    COMPLEMENTARY_6 = "#a98467"
    COMPLEMENTARY_7 = "#6c584c"
    COMPLEMENTARY_8 = "#809bce"
    COMPLEMENTARY_9 = "#95b8d1"
    COMPLEMENTARY_10 = "#b8e0d4"
    COMPLEMENTARY_11 = "#d6eadf"
    COMPLEMENTARY_12 = "#eac4d5"
    LIGHT_1 = "#d6e6ff"
    LIGHT_2 = "#d7f9f8"
    LIGHT_3 = "#ffffea"
    LIGHT_4 = "#fff0d4"
    LIGHT_5 = "#fbe0e0"
    LIGHT_6 = "#e5d4ef"

    @classmethod
    def get_color_list(cls) -> List[str]:
        return [color for color in cls.__members__.keys()]

    @classmethod
    def get_color_codes(cls) -> List[str]:
        return [cls[color].value for color in cls.__members__.keys()]
