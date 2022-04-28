# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 08:15:56 2018

@author: LiWenHuang
"""

from selenium import webdriver
from time import sleep
url = 'ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/XML/'
browser = webdriver.Chrome()
browser.maximize_window()
options = webdriver.ChromeOptions() 
#prefs = {'profile.default_content_settings.popups': 0, 'download.default_directory': 'D:\pubchem'}
#options.add_experimental_option('prefs', prefs) 
#driver = webdriver.Chrome(executable_path='D:\\chromedriver.exe', chrome_options=options) 
browser.get(url)
sleep(3)
for i in range(1,25001,25000): # need to be modified based on the given pubchem xml files
    j = i+24999
    browser.find_element_by_link_text("Compound_"+str("%09d" % i )+"_"+str("%09d" % j )+".xml.gz").click()
    sleep(8)
print("finish")