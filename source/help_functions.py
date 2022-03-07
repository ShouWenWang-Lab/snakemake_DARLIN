import os

import numpy as np
import pandas as pd
from matplotlib import cbook, cm, colors, rcParams


def update_CARLIN_dir(CARLIN_root_folder, template):
    if template == "cCARLIN":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/cCARLIN'"
        )
        # os.system(f"cp {CARLIN_root_folder}/cCARLIN/@CARLIN_def/CARLIN_def_cCARLIN.m {CARLIN_root_folder}/cCARLIN/@CARLIN_def/CARLIN_def.m")
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/cCARLIN"
    elif template == "Tigre":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Tigre_CARLIN'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Tigre_CARLIN"
    elif template == "Tigre_2022":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Tigre_CARLIN_2022'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Tigre_CARLIN_2022"
    elif template == "Rosa":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Rosa_CARLIN'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Rosa_CARLIN"
    else:
        raise ValueError(
            "The input template should be among {Rosa, Tigre_2022, Tigre, cCARLIN}"
        )
    return Actual_CARLIN_dir


def training_notification(
    to_addr="shouwennotification@gmail.com", msg="pythia training finished"
):
    """
    Note that the msg should not contain ':'.
    """
    import smtplib
    import ssl

    smtp_server = "smtp.gmail.com"
    port = 587  # For starttls
    sender_email = "shouwennotification@gmail.com"
    # to_addr='wangsw09@gmail.com'
    password = "as$3+1245QWcn"  # input("Type your password and press enter: ")

    # msg='pythia training finished---dadadg-=++adgag'

    # Create a secure SSL context
    context = ssl.create_default_context()
    print(msg)

    # Try to log in to server and send email
    try:
        server = smtplib.SMTP(smtp_server, port)
        server.ehlo()  # Can be omitted
        server.starttls(context=context)  # Secure the connection
        server.ehlo()  # Can be omitted
        server.login(sender_email, password)
        server.sendmail(sender_email, to_addr, msg)
        # TODO: Send email here
    except Exception as e:
        # Print any error messages to stdout
        print(e)
    finally:
        server.quit()


def set_rcParams(fontsize=12, color_map=None, frameon=None):
    """Set matplotlib.rcParams to cospar defaults."""
    # check here if you want to customize it: https://matplotlib.org/stable/tutorials/introductory/customizing.html

    # dpi options (mpl default: 100, 100)
    rcParams["figure.dpi"] = 100
    rcParams["savefig.dpi"] = 150

    # figure (mpl default: 0.125, 0.96, 0.15, 0.91)
    rcParams["figure.figsize"] = (6, 4)
    # rcParams["figure.subplot.left"] = 0.18
    # rcParams["figure.subplot.right"] = 0.96
    # rcParams["figure.subplot.bottom"] = 0.15
    # rcParams["figure.subplot.top"] = 0.91

    # lines (defaults:  1.5, 6, 1)
    rcParams["lines.linewidth"] = 1.5  # the line width of the frame
    rcParams["lines.markersize"] = 6
    rcParams["lines.markeredgewidth"] = 1

    # font
    rcParams["font.sans-serif"] = [
        "Arial",
        "Helvetica",
        "DejaVu Sans",
        "Bitstream Vera Sans",
        "sans-serif",
    ]

    fontsize = fontsize
    labelsize = 0.92 * fontsize

    # fonsizes (mpl default: 10, medium, large, medium)
    rcParams["font.size"] = fontsize
    rcParams["legend.fontsize"] = labelsize
    rcParams["axes.titlesize"] = fontsize
    rcParams["axes.labelsize"] = labelsize

    # legend (mpl default: 1, 1, 2, 0.8)
    rcParams["legend.numpoints"] = 1
    rcParams["legend.scatterpoints"] = 1
    rcParams[
        "legend.handlelength"
    ] = 1.5  # change it from 1 to 1.5 to allow seaborn function properly
    rcParams["legend.handletextpad"] = 0.4
    rcParams["pdf.fonttype"] = 42

    # color cycle
    # rcParams["axes.prop_cycle"] = cycler(color=vega_10)

    # axes
    rcParams["axes.linewidth"] = 0.8
    rcParams["axes.edgecolor"] = "black"
    rcParams["axes.facecolor"] = "white"

    # ticks (mpl default: k, k, medium, medium)
    rcParams["xtick.color"] = "k"
    rcParams["ytick.color"] = "k"
    rcParams["xtick.labelsize"] = labelsize
    rcParams["ytick.labelsize"] = labelsize

    # axes grid (mpl default: False, #b0b0b0)
    rcParams["axes.grid"] = False
    rcParams["grid.color"] = ".8"

    # color map
    rcParams["image.cmap"] = "Reds" if color_map is None else color_map

    # spines
    rcParams["axes.spines.right"] = False
    rcParams["axes.spines.top"] = False

    # frame (mpl default: True)
    frameon = False if frameon is None else frameon
    global _frameon
    _frameon = frameon
