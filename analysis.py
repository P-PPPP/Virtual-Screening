import gzip, os , requests, json
from gemmi import cif

def Carbon_Sites(uniprot_id:str)-> list:
    # 糖基化位点需要去uniprot找 https://www.uniprot.org/help/api
    # https://rest.uniprot.org/uniprotkb/P05046.txt 中的FT   CARBOHYD        107 代表着107号位置上为糖基化位点
    carbon_sites = []
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt")
    for i in r.text.split("\n"):
        if "FT   CARBOHYD" in i: carbon_sites.append(i.replace("FT   CARBOHYD","").replace(" ","")) #是糖基化位点
    return carbon_sites

def ugzip():
    def un_gz(file_name):
        """ungz zip file"""
        f_name = file_name.replace(".gz", "")
        #获取文件的名称，去掉
        g_file = gzip.GzipFile(file_name)
        #创建gzip对象
        open(f_name, "wb+").write(g_file.read())
        #gzip对象用read()打开后，写入open()建立的文件里。
        g_file.close() #关闭gzip对象
    for i in [fn for fn in os.listdir("./data/") if fn.endswith(".gz")]:
        un_gz(f"./data/{i}")

def readfiles(p:os.PathLike):
    # read and parse a CIF file
    main = []
    for i in [fn for fn in os.listdir(p) if fn.endswith(".cif")]:
        print(i)
        k = {
            "id":i.replace(".cif",""),
        }
        doc = cif.read_string(open(os.path.join(p,i),"r",encoding="utf-8").read())
        block = doc.sole_block()
        uniprot = block.find_value('_struct_ref.pdbx_db_accession')
        
        seq,combine = [],[]
        if block.find_value('_entity_poly.pdbx_seq_one_letter_code') == None:
            # 说明是多链组合的polymer,这时需要loop一下
            _loop = block.find_loop('_entity_poly.pdbx_seq_one_letter_code')
            for i in range(0,len(_loop)): seq.append(_loop[i].replace(";","").replace(" ","").replace("\n",""))
        else: seq = [block.find_value('_entity_poly.pdbx_seq_one_letter_code').replace(";","").replace(" ","").replace("\n","")]

        # for beta sheet
        start = block.find_loop('_struct_sheet_range.beg_label_seq_id') # amino acid begin
        end = block.find_loop('_struct_sheet_range.end_label_seq_id') # amino acid end
        assert len(start) == len(end)
        for i in range(0,len(start)): combine.append([start[i],end[i]]) # seq -1

        try:
            cs = Carbon_Sites(uniprot)
        except: cs = []
        
        k.update({
            "seq":seq,
            "uniprot":uniprot,
            "sheet":combine,
            "Carbon_Sites":cs
        })
        main.append(k)
    open("./data/main.json","w",encoding="utf-8").write(json.dumps(main,indent=2,ensure_ascii=False))

#ugzip()
readfiles("./data")