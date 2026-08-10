from slophep.Core.Parameter import ParameterUser
from slophep.Core.user_registry import FFregistry

import logging
import slophep
slophep.logger.setLevel(logging.ERROR)

transition_dict = {
    "BdToDst"   : r"$B^0 \to D^{*-}$"          ,
    "BdToD"     : r"$B^0 \to D^{-}$"           ,
    "BdToPi"    : r"$B^0 \to \pi^{-}$"         ,
    "BuToD"     : r"$B^\pm \to D^{0}$"         ,
    "BcToJpsi"  : r"$B_c \to J/\psi$"          ,
    "LbToLc"    : r"$\Lambda_b \to \Lambda_c$" ,
    "BdToD0st"  : r"$B^0 \to D_0^{*+}$"        ,
    "BuToD0st"  : r"$B^\pm \to D_0^{*0}$"      ,
    "BdToD1"    : r"$B^0 \to D_1$"             ,
    "BuToD1"    : r"$B^\pm \to D_1^0$"         ,
    "BdToD1st"  : r"$B^0 \to D_1^\prime$"      ,
    "BuToD1st"  : r"$B^\pm \to D_1^{\prime 0}$",
    "BdToD2st"  : r"$B^0 \to D_2^{*}$"         ,
    "BuToD2st"  : r"$B^\pm \to D_2^{*0}$"      ,
    "BsToDs1"   : r"$B^0_s \to D_{s1}$"        ,
    "BsToDs2st" : r"$B^0_s \to D_{s2}$"        ,
    "BsToDsst"  : r"$B^0_s \to D_s^{*-}$"      ,
    "BsToK"     : r"$B^0_s \to K^{-}$"         ,
}


def make_md_table(ff: ParameterUser, toprint: bool):
    table = [
        "| Qualified Name | Value |",
        "|----------------|-------|"
    ]
    if ff is None:
        return []
    params = ff.user_params_defaults()
    for ipar, ival in params.items():
        table.append(f"| `{ipar}` | {ival} |")

    if toprint:
        for iline in table:
            print(iline)

    return table


def print_transition_info(transition: str, transition_label: str):
    print("## "+transition_label)
    schemes = FFregistry.find_all_containing(f"FF{transition}@")
    for ischeme, iff in schemes.items():
        print("### "+ischeme.split("@")[-1])
        table = make_md_table(iff(), True)
        print("\n")


if __name__ == "__main__":
    for itransition, ilabel in transition_dict.items():
        print_transition_info(itransition, ilabel)
