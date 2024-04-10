from selenium import webdriver
from selenium.webdriver.chrome.options import Options
chrome_options = Options()
chrome_options.add_experimental_option("detach", True)
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
# Create a new instance of the Chrome driver
driver = webdriver.Chrome( options=chrome_options)
driver.get("https://www.genome.jp/kegg/mapper/reconstruct.html")

file_input=driver.find_element(By.XPATH,'//input[@type="file"][@name="color_list"]')
file_input.send_keys("/home/davide/Downloads/user_ko_filtered.txt")
submit_button = driver.find_element(By.XPATH, '//input[@type="submit"][@value="Exec"]')
submit_button.click()
# Close the browser window
page_source = driver.page_source

# Write the page source to a file
import os
cwd=os.getcwd()
with open(os.path.join(cwd,"KEGG.html"), "w") as f:
    f.write(page_source)
brite_link = driver.find_element(By.XPATH, '//li[@class="off"]/a[starts-with(@href, "javascript:submit_mapper(\'find_brite_object\',")]')
brite_link.click()
page_source = driver.page_source
