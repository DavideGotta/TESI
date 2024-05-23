import gspread
gc = gspread.oauth()

sh = gc.open("CCMEE_29")
worksheet = sh.sheet1
list_of_lists = worksheet.get_all_values()
worksheet = sh.get_worksheet(1)
#worksheet = sh.add_worksheet(title="CCMEE_29conevidenza", rows="1000", cols="20")
worksheet.update(list_of_lists)
#now highlight in bold first row and first column
worksheet.format("A1:Z1", {"textFormat": {"bold": True}})
worksheet.format("A1:A1000", {"textFormat": {"bold": True}})