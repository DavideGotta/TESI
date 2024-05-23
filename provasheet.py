import os.path
import pandas as pd
from google.auth.transport.requests import Request
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
import math

# If modifying these scopes, delete the file token.json.
SCOPES = ["https://www.googleapis.com/auth/spreadsheets"]

# The ID and range of a sample spreadsheet.
SPREADSHEET_ID = "1jhywU5V8JUCcDWYee7gQlt1h6qSZw9apaFylSdWMqHQ"



def write_to_sheet(service, spreadsheet_id, range_, values):
  request = service.spreadsheets().values().update(
    spreadsheetId=spreadsheet_id,
    range=range_,
    valueInputOption="USER_ENTERED",
    body={"values": values},
  )
  response = request.execute()
  return response

def main():
  """Shows basic usage of the Sheets API.
  Prints values from a sample spreadsheet.
  """
  creds = None
  # The file token.json stores the user's access and refresh tokens, and is
  # created automatically when the authorization flow completes for the first
  # time.
  if os.path.exists("token.json"):
    creds = Credentials.from_authorized_user_file("token.json", SCOPES)
  # If there are no (valid) credentials available, let the user log in.
  if not creds or not creds.valid:
    if creds and creds.expired and creds.refresh_token:
      creds.refresh(Request())
    else:
      flow = InstalledAppFlow.from_client_secrets_file(
          "credentials.json", SCOPES
      )
      creds = flow.run_local_server(port=0)
    # Save the credentials for the next run
    with open("token.json", "w") as token:
      token.write(creds.to_json())

  try:
    service = build("sheets", "v4", credentials=creds)

    # Call the Sheets API
    SPREADSHEET_ID = "1Sn7Gr2YHGeP7Y7Lb-WqNVVj5M6xycWZECADkJEG6Tzc"
    a, b, dx = 10, 39, 0.5

    def f(x):
      return math.sin(math.sqrt(x)) / math.sqrt(x)

    values = [["step", "intervallo", "punto medio", "Area rettangolo"]]

    values += [[i+1, f"{a+i*dx}-{a+(i+1)*dx}", a+i*dx+dx/2, f(a+i*dx+dx/2)*dx] for i in range(int((b - a) / dx))]

    values.append(["", "", "Risultato", sum(row[3] for row in values[1:])])

    write_to_sheet(service, SPREADSHEET_ID, "Foglio1!A1:Z100", values)


    #print(values)
    #create a dataframe from the values

  except HttpError as err:
    print(err)


if __name__ == "__main__":
  main()



