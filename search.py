import json, itertools 

conditions = {
    'sequence_part' : ["SHD","HHH"],
    'Structure': ["sheet"]
}

main = json.loads(open("./data/main.json","r",encoding="utf-8").read())

def search_main(item:dict)->dict:
    result = {
        "result":False,
        "info":[],
        "id":item["id"]
    }
    for z,k in itertools.product(conditions["sequence_part"],item["seq"]):
        # compare sequences
        info = {
            "positions":[],
            "betaSheet":[]
        }
        if k.find(z) != -1:
            result["result"] = True
            # 目标短肽在序列中
            info.update({"seq":k,})
            info["positions"] =  [i+1 for i in range(len(k)) if k[i:i+len(z)] == z]
        if len(info["positions"]) != 0:
            for a,b in itertools.product(info["positions"],item["sheet"]):
                if  a - int(b[0]) > 0 and int(b[1]) -a > 0: #在β折叠区域
                    info["betaSheet"].append(a)
        result["info"].append(info)
    return result

if __name__ == "__main__":
    q = []
    for i in main:
        result = search_main(i)
        if result["result"] == True: q.append(result)

    open("./result.json","w",encoding="utf-8").write(json.dumps(q,indent=1,ensure_ascii=False))