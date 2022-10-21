library(shinyobjects)

# widget view
pickerInput("cellInput","Cells",
            choices=c("Acinar","Alpha","Beta"),
            options = list(`actions-box` = TRUE), multiple = T
) %>% view_ui()


tbl_out # chr, start,end, proxy, ocr_chr, ocr_start, ocr_end, id, gene_name, gene_id cell, called, geneOCR_table_type


##################### test functions #####################
### read_snipa
file="test/proxySearch.results.csv"
proxy_df = read_snipa(file)

### cal_proxy
proxy_search = T
sentinel = "rs114526150"
window_size = 5e5
pop = "EUR"
r2 = 0.8
token = "9870c4774bad"

proxy_df2 = cal_proxy(proxy_search, sentinel, window_size, pop, r2, token)

cal_proxy(F, sentinel, window_size, pop, r2, token)

### v2g_mapping_batch
bi_df = v2g_mapping_batch(proxy_df=proxy_df2[["proxies"]], geneOCR_table_type="Binary", cells=c("Acinar","Alpha","Beta"))
ABC_df = v2g_mapping_batch(proxy_df=proxy_df2[["proxies"]], geneOCR_table_type="ABC model", cells=c("Acinar","Alpha","Beta"))
promoter_df = v2g_mapping_batch(proxy_df=proxy_df2[["proxies"]], geneOCR_table_type="Promoter", cells=c("Acinar","Alpha","Beta"))
input = list(
        hicInput=c("Binary", "ABC model", "Promoter"),
        cellInput=c("Acinar","Alpha","Beta")
)
tbl_out =  do.call(
        "bind_rows",
        lapply(
                input$hicInput, # geneOCR_table_type
                v2g_mapping_batch,
                proxy_df = proxy_df2[["proxies"]],
                cells = input$cellInput
        )
)

### plot_snp_gene
plot_snp_gene(tbl_out)
plot_snp_gene(tbl_out %>% filter(geneOCR_table_type="ABC model"))
plot_snp_gene(tbl_out %>% mutate(geneOCR_table_type="Promoter"))

### search_geneOC
gene_id = gene_type %>%
        filter(gene_name=="IGF2") %>%
        distinct(gene_id) %>%
        dplyr::slice(1) %>%
        pull(gene_id)
search_geneOCR(
        gene_id,
        geneOCR_table_type="Binary",
        cells=c("Acinar","Alpha","Beta")
)

gene_id = gene_type %>%
        filter(gene_name=="FAM87B") %>%
        distinct(gene_id) %>%
        dplyr::slice(1) %>%
        pull(gene_id)
search_geneOCR(
        gene_id,
        geneOCR_table_type="ABC model",
        cells=c("Acinar","Alpha","Beta")
)

gene_id = gene_type %>%
        filter(gene_name=="IGF2") %>%
        distinct(gene_id) %>%
        dplyr::slice(1) %>%
        pull(gene_id)
search_geneOCR(
        gene_id,
        geneOCR_table_type="Promoter",
        cells=c("Acinar","Alpha","Beta")
)


input = list(
        hicInput=c("Binary", "ABC model", "Promoter"),
        cellInput=c("Acinar","Alpha","Beta")
)
tbl_out =  do.call(
        "bind_rows",
        lapply(
                input$hicInput, # geneOCR_table_type
                search_geneOCR,
                gene_id=gene_id,
                cells = input$cellInput
        )
)
plot_ocr_gene(tbl_out, plot_type="arc")
####################### input ###########################
dummy_input <- list(
        button = "",
        gene_name = "",
        hicInput = "",
        pop = "",
        precal = "",
        proxy_file = "",
        proxy_search = "",
        r2 = "",
        searchItem = "",
        sentinel = "",
        token = "",
        window_size = ""
)

###################### additional server  code ##################3
# apply gene filter
input = list(
        gene_filter=T,
        open_gene_only=T

 )
# input$gene_filter (T)
# input$open_gene_only (T)

if (input$gene_filter==T){
        tbl_out = tbl_out %>%
                semi_join(
                        gene_type %>%
                                filter(gene_type=="protein_coding")
                )
}

if (input$open_gene_only==T){
        tbl_out = tbl_out %>%
                semi_join(
                        refOCR_cell_simple %>%
                                filter(!is.na(gene_ids)) %>%
                                separate_rows(gene_ids, sep="\\|") %>%
                                dplyr::rename(gene_id=gene_ids)
                )
}

        
