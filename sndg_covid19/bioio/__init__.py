def country_from_gisaid(fasta_id):
    rid = fasta_id.replace("hCoV-19/", "")

    code, gisaid, sdate = rid.split("|")
    code_len = len(code.split("/"))
    assert code_len in [2,3, 4], code_len
    if code_len in [2,3]:
        # hCoV-19/Wuhan-Hu-1/2019|EPI_ISL_402125|2019-12-31 can be 2
        country = code.split("/")[0]
    elif code_len == 4:
        # hCoV-19/tiger/USA/NY-040420/2020|EPI_ISL_420293|2020-04-02 can be 4
        country = code.split("/")[1]
    return country
