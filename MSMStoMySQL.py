# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 16:51:05 2018

@author: LiWenHuang
"""

import xml.etree.ElementTree as ET
##連接mysql、建立表格與欄位
import mysql.connector as mydb 
conn = mydb.connect(host = '127.0.0.1',user = 'root',passwd = '1234',
    db = 'compound',port = 3306,charset = 'utf8')
cursor = conn.cursor()
print("Connect Mysql")
tables = '''create table msmsdata20181219(collection_date text ,collision_energy_level text ,
                                      collision_energy_voltage text ,created_at text not null,
                                      database_id text not null,id text not null,inchi_key text ,
                                      instrument_type text not null,ionization_mode text not null,
                                      molecule_id text ,peak_counter text not null,sample_assessment text ,sample_concentration text,
                                      sample_concentration_units text ,sample_mass text ,ref_text text ,
                                      sample_mass_units text ,spectra_assessment text ,updated_at text not null,
                                      database_source text not null,database_id0 text  ,id0 text ,pubmed_id text,
                                      spectra_id text not null,spectra_type text not null,ms2_peak MEDIUMTEXT not null)'''
cursor.execute(tables)
print("Table Set OK")
def read_compound_data(filename):
    data = []
    subdata = []
    peakdata = []
    combine_subdata=[]
    label1 = ['collection-date','collision-energy-level','collision-energy-voltage','created-at','database-id',
             'id','instrument-type','ionization-mode','peak-counter','sample-assessment','sample-concentration',
             'sample-concentration-units','sample-mass','sample-mass-units','spectra-assessment','updated-at']
    label2 = ['ref-text','database','database-id','id','pubmed-id','spectra-id','spectra-type']
    tree = ET.parse(filename)
    root = tree.getroot()        
    for i in range(0,len(label1),1):
        if len(root.findall(label1[i])) == 0 and root.text.strip()=='' :
                text="no exist"
                data.append(text)
        else:
                for subroot in root.findall(label1[i]):
                    if len(subroot)==0 and subroot.text==None:
                        text='None'
                        data.append(text)
                    elif len(subroot) == 0 and subroot.text.strip()=='':
                        text="NULL"
                        data.append(text)
                    elif len(subroot)==0 and subroot.text.strip():
                        text=subroot.text
                        data.append(text)
    for j in range(0,len(label2),1):
        for REF_L1 in root.findall('references'):
            if len(REF_L1) == 0 :
                text="no exist"
                subdata.append(text)
            else:
                for REF_L2 in REF_L1:
                    for REF_L3 in REF_L2.findall(label2[j]) :
                        if len(REF_L3)==0 and REF_L3.text==None:
                            text='None'
                            subdata.append(text)
                        elif len(REF_L3) == 0 and REF_L3.text.strip()=='':
                            text="NULL"
                            subdata.append(text)
                        elif len(REF_L3)==0 and REF_L3.text.strip():
                            text=REF_L3.text                
                            subdata.append(text)    
    if len(root.findall('ms-ms-peaks')) == 0 and root.text.strip()=='' :
        text="no exist"
        subdata.append(text)
    else:
        for peak_L1 in root.findall('ms-ms-peaks'):
            for peak_L2 in peak_L1:
                d1 = peak_L2[0].text
                d2 = peak_L2[3].text
                d3 = peak_L2[2].text
                d4 = peak_L2[1].text
                peakdata = [d1,d2,d3,d4]
                combine_peakdata = ";".join(peakdata)
                combine_subdata.append(combine_peakdata)
            combine_peakdata2 = ";".join(combine_subdata)
            subdata.append('\"'+combine_peakdata2+'\"')
    data.insert(6,'')
    data.insert(9,'')
    data.insert(15,'\"'+str(subdata[0])+'\"')
    if len(data)+len(subdata) != 27:
        print(filename)
    else:
        for n in range(1,8,1):
            data.insert(n+18,subdata[n])
        try:
            cursor.execute("insert into msmsdata20181219 values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"\
                           ,[data[i] for i in range(len(data))])
    data.clear()
    subdata.clear()
    peakdata.clear()
    combine_subdata.clear()
#開始讀檔
count = 0
for i in range(1,23153,1):
    count += 1
    file = "01 ("+str("%d" % i )+")"+".xml"
    read_compound_data(file)
    if count%2500 == 0: 
        print(count,'ok')
print("database finish")
conn.close()