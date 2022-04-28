# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 19:51:01 2018

@author: LiWenHuang
"""
import os
import os.path
import gzip
import xml.etree.ElementTree as ET
#連接mysql、建立表格與欄位
import mysql.connector as mydb 
conn = mydb.connect(host = '127.0.0.1',user = 'root',passwd = '1234',
    db = 'compound',port = 3306,charset = 'utf8')
cursor = conn.cursor()
print("Connect Mysql")
tables = '''create table compounddata20181119(NO INT not null PRIMARY KEY,IUPACName_Traditional text not null,
                                  Mass_Exact DOUBLE not null,MolecularFormula text not null,SMILES_Canonical text not null,
                                  SMILES_Isomeric text not null,Weight_MonoIsotopic DOUBLE not null)ENGINE = MyISAM'''
cursor.execute(tables)
print("Table Set OK")
#Def1:建立儲存錯誤檔案名稱
def storage_errorfile(path):
    errorfile.append(path)
#開始讀檔與分析資料
#Def2:分析與建立架構
def read_compound_data(filename):
    label = '{http://www.ncbi.nlm.nih.gov}' #xml檔裡所有標籤值前方字串
    tree = ET.parse(filename)
    root = tree.getroot()
    data = {} 
    compounddata = []
    for subroot in root:
        for ID_L1 in subroot.findall(label+'PC-Compound_id'):
            IDtext = ID_L1[0][0][0].text
            compounddata.append(IDtext)
        for Urn_L1 in subroot.findall(label+'PC-Compound_props'):
            for Urn_L2 in Urn_L1:
                a = Urn_L2[0][0][0].text
                b = Urn_L2[0][0][1].text
                c = Urn_L2[1][0].text
                if b.strip() == '':
                    data.setdefault(a,c)
                else:
                    data.setdefault((a,b),c)
            d1 = '\"'+str(data.get(('IUPAC Name', 'Traditional')))+'\"'
            d2 = float(data.get(('Mass', 'Exact')))
            d3 = data.get('Molecular Formula')
            d4 = data.get(('SMILES', 'Canonical'))
            d5 = data.get(('SMILES', 'Isomeric'))
            d6 = float(data.get(('Weight', 'MonoIsotopic')))
            #多執行序
            urndata = [d1,d2,d3,d4,d5,d6]
            compounddata += urndata
        cursor.execute("insert into compounddata20181119 values (%s,%s,%s,%s,%s,%s,%s)",[compounddata[i] for i in range(len(compounddata))])
        urndata.clear()
        compounddata.clear()
        data.clear()
    print(count,',',filename,"ok")
#Def3:讀壓縮檔
def read_gz_file(path):
    if os.path.exists(path):
        try:
            pubchemdata = gzip.open(path, 'r')
            read_compound_data(pubchemdata)
        except EOFError :
            storage_errorfile(path)
            print(path,"Error!FILE IS LOST")
        except FileNotFoundError:
            storage_errorfile(path)
            print("lostconnect D:// ",path)
    else:
        storage_errorfile(path)
        print(count,'the path [{}] is not exist!'.format(path))
#開始讀檔    
errorfile =[]
count = 0
for i in range(1,134825001,25000):
    count += 1
    j = i+24999
    file = "Compound_"+str("%09d" % i )+"_"+str("%09d" % j )+".xml.gz"
    read_gz_file('D:\compound'+'\\'+file)
print(errorfile)
print("database finish")

#['D:\\compound\\Compound_086350001_086375000.xml.gz', 'D:\\compound\\Compound_086475001_086500000.xml.gz']