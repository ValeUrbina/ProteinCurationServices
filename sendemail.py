import sys
import os
import email
import smtplib
import ssl
from email import encoders
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

sender_email = "proteincurationapp@gmail.com"
password = "tesis2proteincurationapppassword"


def main(PATH, operation_id, receiver_mail, filepath, pfam_code):

    receiver_email = receiver_mail

    message = MIMEMultipart("alternative")
    message["Subject"] = "Protein Curation App - Operation ID"
    message["From"] = sender_email
    message["To"] = receiver_email

    # Create the plain-text and HTML version of your message
    text = """\
    Hi,
    Here's your Operation ID: """ + operation_id + """ 
    Don't forget to enter the code to resume work"""
    html = """\
    <html>
    <body>
        <p>Dear curator,
        <br>Here's your Operation ID: """ + operation_id + """<br>
        <br>Don't forget to enter the code to resume work<br>
        </p>
    </body>
    </html>
    """

    # Turn these into plain/html MIMEText objects
    part1 = MIMEText(text, "plain")
    part2 = MIMEText(html, "html")

    # Add HTML/plain-text parts to MIMEMultipart message
    # The email client will try to render the last part first
    message.attach(part1)
    message.attach(part2)

    # Send the seed - Document directory
    filename = os.path.join(PATH, 'pfam_data', filepath, pfam_code, 'SEED')

    # Open SEED file in binary mode
    with open(filename, "rb") as attachment:
        # Add file as application/octet-stream
        # Email client can usually download this automatically as attachment
        part = MIMEBase("application", "octet-stream")
        part.set_payload(attachment.read())

    # Encode file in ASCII characters to send by email
    encoders.encode_base64(part)

    # Add header as key/value pair to attachment part
    part.add_header(
        "Content-Disposition",
        f"attachment; filename= {pfam_code + '-SEED'}",
    )

    # Add attachment to message and convert message to string
    message.attach(part)
    text = message.as_string()

    # Log in to server using secure context and send email
    context = ssl.create_default_context()
    with smtplib.SMTP_SSL("smtp.gmail.com", 465, context=context) as server:
        server.login(sender_email, password)
        server.sendmail(sender_email, receiver_email, text)

    return 1


if __name__ == "__main__":
    main()
