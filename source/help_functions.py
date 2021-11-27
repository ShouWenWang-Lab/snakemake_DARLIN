import os
import pandas as pd
import numpy as np


def update_CARLIN_dir(CARLIN_root_folder,template):
    if template=='cCARLIN':
        os.system(f"rsync -avP {CARLIN_root_folder}/Custom_CARLIN/ {CARLIN_root_folder}/cCARLIN")
        os.system(f"cp {CARLIN_root_folder}/cCARLIN/@CARLIN_def/CARLIN_def_cCARLIN.m {CARLIN_root_folder}/cCARLIN/@CARLIN_def/CARLIN_def.m")
        Actual_CARLIN_dir=f"{CARLIN_root_folder}/cCARLIN"
    else:
        os.system(f"rsync -avP {CARLIN_root_folder}/Custom_CARLIN/ {CARLIN_root_folder}/Tigre_CARLIN")
        os.system(f"cp {CARLIN_root_folder}/Tigre_CARLIN/@CARLIN_def/CARLIN_def_Tigre.m {CARLIN_root_folder}/Tigre_CARLIN/@CARLIN_def/CARLIN_def.m")
        Actual_CARLIN_dir=f"{CARLIN_root_folder}/Tigre_CARLIN"
        
    return Actual_CARLIN_dir


def training_notification(
    to_addr="shouwennotification@gmail.com", msg="pythia training finished"
):
    """
    Note that the msg should not contain ':'.
    """
    import smtplib, ssl

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