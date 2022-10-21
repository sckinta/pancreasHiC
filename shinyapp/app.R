library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinythemes)
library(shinyjs)
library(vroom)
# library(data.table)
# options(tidyverse.quiet = TRUE)
library(dplyr)
library(tidyr)
library(ggplot2)
options(dplyr.summarise.inform = FALSE)
library(myvariant)
library(LDlinkR)
library(gridExtra)
# library(GenomicRanges)
library(DFbedtools)
library(igraph)
library(ggraph)



############################ load database ####################
gene_type = vroom("data/gencodeV19_gene_type.protein_lncRNA.txt.gz")
genenames = gene_type %>% distinct(gene_name) %>% pull(gene_name)
ABCresults_cellHiC_simple = vroom("data/ABCresults_cellHiC_simple.txt.gz")
geneOCR_hicNaive_simple = vroom("data/geneOCR_hicNaive_simple.txt.gz")
refOCR_cell_simple = vroom("data/refOCR_cell_simple2.txt.gz")

########################### functions ##########################
read_snipa <- function(proxy_file){
    message("Read in proxy file in SNiPA output format...<br>")
    proxies = suppressMessages(vroom(proxy_file)) %>%
        dplyr::select(query = QRSID, proxy=RSID, proxy_chr=CHR, proxy_pos=POS2, R.squared=R2) %>%
        mutate(proxy_chr=gsub("^","chr",proxy_chr)) %>%
        mutate(proxy = ifelse(!grepl("rs",proxy), glue::glue("{proxy_chr}:{proxy_pos}"), proxy))
   file=NULL
    if(nrow(proxies) > 0){
        list(proxies=proxies, error=NULL)
    }else{
        list(proxies=NULL, error="no proxies found in given file")
    }
}

cal_proxy <- function(proxy_search, sentinel, window_size, pop, r2, token){
    error = 0
    if (proxy_search){
        # with proxy_search
        message("Calculate proxies...<br>")

        # format input sentinel
        if (!grepl("rs",sentinel)){
            if (grepl(":", sentinel)){
                sentinel = tolower(sentinel)
                if(!grepl("chr",sentinel)){
                    sentinel=gsub("^","chr",sentinel)
                }
            }else{
                error="input sentinel does not match required format; format need to be rsID or chromosome coordinate (e.g. chr7:24966446)"
            }
        }

        # LDproxy
        proxies = LDproxy(
            snp = sentinel, pop = pop,
            r2d = "r2", token = token
        )
        if(grepl("error",proxies[1, 1])){
            error=proxies[1, 1]
        }

        # clean LDproxy output
        if(error==0){
            proxies = proxies %>% as_tibble() %>%
                mutate(query=sentinel) %>%
                mutate_if(function(x){is.factor(x)}, function(x){as.character(x)}) %>%
                dplyr::select(query, proxy=RS_Number, proxy_coord=Coord, R.squared=R2, Distance) %>%
                separate(proxy_coord, c("proxy_chr","proxy_pos"),sep=":") %>%
                mutate(proxy_pos=as.integer(proxy_pos)) %>%
                mutate(proxy = ifelse(!grepl("rs",proxy), glue::glue("{proxy_chr}:{proxy_pos}"), proxy))
            if(nrow(proxies %>% filter(Distance==0))!=1){
                error="proxy search cannot be done for given sentinel; try to map sentinel by un-checking proxy search"
            }else{
                proxies = proxies %>%
                    filter(R.squared >= r2, abs(Distance) <= window_size) %>%
                    dplyr::select(-Distance)
                if (nrow(proxies)==0){
                    error="no proxies passed R2 threshold; try to lower R2 and run again"
                }
            }
        }
    }else{
        # no proxy search
        message("Search input sentinel genomic location...<br>")

        if (grepl("rs",sentinel)){
            # queryVariant{myvariant}
            queryResults = queryVariant(sentinel,  return.as="DataFrame")
            chr = paste0("chr",queryResults$hits$chrom)
            pos = queryResults$hits$hg19$end
        } else if (grepl(":", sentinel)){
            tmp = strsplit(sentinel,":")[[1]]
            chr = gsub("chr|Chr|CHR|CHROM","",tmp[1])
            chr = gsub("^","chr",chr)
            pos = as.integer(tmp[2])
        }
        if (is.null(pos) | chr=="chr"){
            error="cannot find position for given snp in database"
        }else{
            proxies = tibble(
                query=sentinel,
                proxy=sentinel,
                proxy_chr=chr,
                proxy_pos=pos,
                R.squared=1
            ) %>% distinct()
        }
    }

    # return output
    if (error==0){
        list(proxies=proxies, error=NULL)
    }else{
        list(proxies=NULL, error=error)
    }
}

v2g_mapping_batch <- function(proxy_df, geneOCR_table_type, cells){
    # proxy_df format
    proxy_df = proxy_df %>%
        mutate(start=proxy_pos-1) %>%
        dplyr::select(chr=proxy_chr, start, end=proxy_pos, proxy) %>%
        distinct()

    if (geneOCR_table_type=="Binary"){
        geneOCR_table = geneOCR_hicNaive_simple %>%
            filter(cell %in% cells)

        # overlap with OCR
        proxy2geneOCR = overlap_df(
            proxy_df,
            geneOCR_table,
            df1_0base=T, df2_0base=T, minoverlap=1L
        )
        proxy2geneOCR = bind_cols(
            proxy2geneOCR$overlap_df1,
            proxy2geneOCR$overlap_df2 %>%
                dplyr::rename(
                    ocr_chr=chr,
                    ocr_start=start,
                    ocr_end=end
                )
        )

        if (nrow(proxy2geneOCR)==0){
            proxy2geneOCR = tibble(
                chr = character(),
                start = integer(),
                end = integer(),
                proxy = character(),
                ocr_chr = character(),
                ocr_start = integer(),
                ocr_end = integer(),
                id = integer(),
                gene_name = character(),
                gene_id = character(),
                cell = character(),
                hicLoops = character(),
                called = integer(),
                geneOCR_table_type = character()
            )
        }else{
            # overlap with OCREnd HiC
            proxy2ocrEnd = overlap_df(
                proxy_df,
                proxy2geneOCR %>%
                    dplyr::select(ocrEnd, resol) %>%
                    distinct() %>%
                    separate(ocrEnd, c("ocrEnd_chr","ocrEnd_start","ocrEnd_end"), remove=F) %>%
                    mutate_at(c("ocrEnd_start","ocrEnd_end"), as.integer) %>%
                    dplyr::select(ocrEnd_chr, ocrEnd_start, ocrEnd_end, ocrEnd, resol),
                df2_chr_col="ocrEnd_chr", df2_start_col="ocrEnd_start", df2_end_col="ocrEnd_end", df1_0base=T, df2_0base=T, minoverlap=1L
            )
            proxy2ocrEnd = bind_cols(
                proxy2ocrEnd$overlap_df1,
                proxy2ocrEnd$overlap_df2 %>%
                    dplyr::select(ocrEnd,resol)
            )
            proxy2geneOCR = semi_join(
                proxy2geneOCR,
                proxy2ocrEnd
            )
            proxy2geneOCR = proxy2geneOCR %>%
                mutate(hicLoop=glue::glue("{resol}+{ocrEnd}+{geneEnd}")) %>%
                dplyr::select(-resol, -ocrEnd, -geneEnd) %>%
                group_by(chr, start, end, proxy, ocr_chr, ocr_start, ocr_end, id, gene_name, gene_id, cell) %>%
                summarise(
                    hicLoops=paste(hicLoop, collapse="|"),
                    called=sum(called)
                ) %>%
                ungroup() %>%
                mutate(called=ifelse(called >= 1, 1, 0)) %>%
                mutate(geneOCR_table_type="Binary")
        }
    }else if (geneOCR_table_type=="ABC model"){

        geneOCR_table = ABCresults_cellHiC_simple %>%
            filter(cell %in% cells)

        # overlap with OCR
        proxy2geneOCR = overlap_df(
            proxy_df,
            geneOCR_table,
            df1_0base=T, df2_0base=T, minoverlap=1L
        )
        proxy2geneOCR = bind_cols(
            proxy2geneOCR$overlap_df1,
            proxy2geneOCR$overlap_df2 %>%
                dplyr::rename(
                    ocr_chr=chr,
                    ocr_start=start,
                    ocr_end=end
                )
        )

        if (nrow(proxy2geneOCR)==0){
            proxy2geneOCR = tibble(
                chr = character(),
                start = integer(),
                end = integer(),
                proxy = character(),
                ocr_chr = character(),
                ocr_start = integer(),
                ocr_end = integer(),
                id = integer(),
                gene_name = character(),
                gene_id = character(),
                cell = character(),
                hicLoops = character(),
                called = integer(),
                geneOCR_table_type = character()
            )
        }else{
            proxy2geneOCR = proxy2geneOCR %>%
                dplyr::rename(called=ABC.Score) %>%
                mutate(geneOCR_table_type="ABC model", hicLoops=NA)

        }
    }else if (geneOCR_table_type=="Promoter"){
        geneOCR_table = refOCR_cell_simple %>%
            filter(cell %in% cells) %>%
            filter(!is.na(gene_ids))
        proxy2geneOCR = overlap_df(
            proxy_df,
            geneOCR_table,
            df1_0base=T, df2_0base=T, minoverlap=1L
        )

        proxy2geneOCR = bind_cols(
            proxy2geneOCR$overlap_df1,
            proxy2geneOCR$overlap_df2 %>%
                dplyr::rename(
                    ocr_chr=chr,
                    ocr_start=start,
                    ocr_end=end
                )
        )

        if (nrow(proxy2geneOCR)==0){
            proxy2geneOCR = tibble(
                chr = character(),
                start = integer(),
                end = integer(),
                proxy = character(),
                ocr_chr = character(),
                ocr_start = integer(),
                ocr_end = integer(),
                id = integer(),
                gene_name = character(),
                gene_id = character(),
                cell = character(),
                hicLoops = character(),
                called = integer(),
                geneOCR_table_type = character()
            )
        }else{
            proxy2geneOCR = proxy2geneOCR %>%
                dplyr::select(chr, start,end, proxy, ocr_chr, ocr_start, ocr_end, id, gene_ids, cell, called) %>%
                separate_rows(gene_ids,sep="\\|") %>%
                dplyr::rename(gene_id=gene_ids) %>%
                left_join(gene_type %>% distinct(gene_id, gene_name)) %>%
                mutate(geneOCR_table_type="Promoter", hicLoops=NA)
        }
    }

    proxy2geneOCR # chr, start,end, proxy, ocr_chr, ocr_start, ocr_end, id, gene_name, gene_id cell, called, geneOCR_table_type
}

search_geneOCR <- function(gene_id, geneOCR_table_type, cells){
    gene_df = gene_type %>% filter(gene_id==!!gene_id) %>%
        dplyr::select(chr=pro_chr, start=pro_start, end=pro_end, gene_id, gene_name)

    if (geneOCR_table_type=="Binary"){
        geneOCR_table = geneOCR_hicNaive_simple %>%
            semi_join(
                gene_df %>%
                    dplyr::select(gene_id, gene_name)
            ) %>%
            filter(cell %in% cells)

        if (nrow(geneOCR_table)==0){
            geneOCR_out = tibble(
                chr=character(),
                start=integer(),
                end=integer(),
                id = integer(),
                cell = character(),
                called = integer(),
                gene_id = character(),
                gene_name = character(),
                geneOCR_table_type = character(),
                hicLoops = character()
            )
        }else{
            geneOCR_out = geneOCR_table %>%
                mutate(hicLoop=glue::glue("{resol}+{ocrEnd}+{geneEnd}")) %>%
                dplyr::select(-resol, -ocrEnd, -geneEnd) %>%
                distinct(chr, start, end, id, cell, called, hicLoop) %>%
                group_by(chr, start, end, id, cell) %>%
                summarise(
                    hicLoops=paste(hicLoop, collapse="|"),
                    called=sum(called)
                ) %>%
                ungroup() %>%
                mutate(called=ifelse(called > 0, 1, 0)) %>%
                bind_rows(
                    gene_df %>%
                        dplyr::slice(rep(1:n(), each=length(cells)))
                ) %>%
                mutate(geneOCR_table_type="Binary")
        }
    }else if (geneOCR_table_type=="ABC model"){
        geneOCR_table = ABCresults_cellHiC_simple %>%
            semi_join(
                gene_df %>%
                    dplyr::select(gene_id, gene_name)
            ) %>%
            filter(cell %in% cells)

        if (nrow(geneOCR_table)==0){
            geneOCR_out = tibble(
                chr=character(),
                start=integer(),
                end=integer(),
                id = integer(),
                cell = character(),
                called = numeric(),
                gene_id = character(),
                gene_name = character(),
                geneOCR_table_type = character(),
                hicLoops = character()
            )
        }else{
            geneOCR_out = geneOCR_table %>%
                distinct(chr, start, end, id, cell, ABC.Score) %>%
                group_by(chr, start, end, id, cell) %>%
                dplyr::slice(which.max(ABC.Score)) %>%
                ungroup() %>%
                dplyr::rename(called=ABC.Score) %>%
                bind_rows(
                    gene_df %>%
                        dplyr::slice(rep(1:n(), each=length(cells)))
                ) %>%
                mutate(geneOCR_table_type="ABC model", hicLoops=NA)
        }
    }else if (geneOCR_table_type=="Promoter"){
        geneOCR_table = refOCR_cell_simple %>%
            filter(!is.na(gene_ids)) %>%
            separate_rows(gene_ids, sep="\\|") %>%
            dplyr::rename(gene_id=gene_ids) %>%
            semi_join(
                gene_df %>%
                    dplyr::select(gene_id, gene_name)
            ) %>%
            filter(cell %in% cells)

        if (nrow(geneOCR_table)==0){
            geneOCR_out = tibble(
                chr=character(),
                start=integer(),
                end=integer(),
                id = integer(),
                cell = character(),
                called = numeric(),
                gene_id = character(),
                gene_name = character(),
                geneOCR_table_type = character(),
                hicLoops = character()
            )
        }else{
            geneOCR_out = geneOCR_table %>%
                distinct(chr, start, end, id, cell, called) %>%
                group_by(chr, start, end, id, cell) %>%
                summarise(called=sum(called)) %>%
                ungroup() %>%
                mutate(called=ifelse(called > 0, 1, 0)) %>%
                bind_rows(
                    gene_df %>%
                        dplyr::slice(rep(1:n(), each=length(cells)))
                ) %>%
                mutate(geneOCR_table_type="Promoter", hicLoops=NA)
        }
    }
    geneOCR_out
}

plot_snp_gene <- function(tbl_out){
    df = tbl_out %>%
        distinct(cell, gene_id, gene_name, called, geneOCR_table_type)
    if (length(unique(df$geneOCR_table_type)) == 1 && unique(df$geneOCR_table_type)!="Binary"){
        if (unique(df$geneOCR_table_type)=="ABC model"){
            df = df %>%
                filter(geneOCR_table_type %in% c("ABC model")) %>%
                group_by(gene_name, gene_id, cell) %>%
                dplyr::slice(which.max(called)) %>%
                ungroup()

            gene_n=length(unlist(df %>% distinct(gene_name)))
            cell_n=length(unlist(df %>% distinct(cell)))

            if (gene_n >= 20){
                df = df %>%
                    semi_join(
                        df %>%
                            group_by(gene_id) %>%
                            summarise(sd=sd(called)) %>%
                            ungroup() %>%
                            arrange(desc(sd)) %>%
                            head(n=20)
                    )
            }

            p = df %>%
                mutate(called=ifelse(called < 0.02, NA, called)) %>%
                ggplot(aes(x=gene_name, y=cell, fill=called)) +
                geom_tile(color="grey50", size=1) +
                scale_fill_gradient(low="skyblue",high="skyblue4",na.value ="white") +
                theme(
                    panel.background = element_blank(),
                    axis.title.x=element_blank(),
                    axis.ticks = element_blank(),
                    axis.text.x = element_text(angle = 90, color = "grey50", size=7),
                    axis.title.y=element_blank()
                ) +
                labs(fill="ABC score")

            if (gene_n >= 20){
                p = p + ggtitle("Most cell variable Gene Top 20")
                gene_n=length(unlist(df %>% distinct(gene_name)))
            }

        }else if (unique(df$geneOCR_table_type)=="Promoter"){
            df = tbl_out %>%
                # filter(geneOCR_table_type=="Promoter") %>%
                left_join(
                    refOCR_cell_simple %>%
                        dplyr::select(id, cell, activity_base)
                ) %>%
                distinct(cell, gene_id, gene_name, called, activity_base)

            gene_n=length(unlist(df %>% distinct(gene_name)))
            cell_n=length(unlist(df %>% distinct(cell)))

            if (gene_n >= 20){
                df = df %>%
                    semi_join(
                        df %>%
                            group_by(gene_id) %>%
                            summarise(sd=sd(called)) %>%
                            ungroup() %>%
                            arrange(desc(sd)) %>%
                            head(n=20)
                    )
            }

            p = df %>%
                mutate(activity_base = ifelse(called==0, NA, activity_base)) %>%
                ggplot(aes(x=gene_name, y=cell, fill=activity_base)) +
                geom_tile(color="grey50", size=1) +
                scale_fill_gradient(low="skyblue",high="skyblue4",na.value ="white") +
                theme(
                    panel.background = element_blank(),
                    axis.title.x=element_blank(),
                    axis.ticks = element_blank(),
                    axis.text.x = element_text(angle = 90, color = "grey50", size=7),
                    axis.title.y=element_blank()
                ) +
                labs(fill="OCR enhancer activity")

            if (gene_n >= 20){
                p = p + ggtitle("Most cell variable Gene Top 20")
                gene_n=length(unlist(df %>% distinct(gene_name)))
            }
        }
    }else{
        df = df %>%
            mutate(called=ifelse(called >= 0.02, 1, 0))
        df = df %>%
                group_by(gene_name, gene_id, cell) %>%
                summarise(
                    called=sum(called),
                    geneOCR_table_type=paste(unique(geneOCR_table_type), collapse = ",")
                ) %>%
                ungroup() %>%
                mutate(called = ifelse(called != 0, 1, 0))

        df = df %>%
            mutate(geneOCR_table_type = ifelse(called==0, NA, geneOCR_table_type))

        gene_n=length(unlist(df %>% distinct(gene_name)))
        cell_n=length(unlist(df %>% distinct(cell)))

        if (gene_n >= 20){
            df = df %>%
                semi_join(
                    df %>%
                        group_by(gene_id) %>%
                        summarise(sd=sd(called)) %>%
                        ungroup() %>%
                        arrange(desc(sd)) %>%
                        head(n=20)
                )
        }

        p = df %>%
            ggplot(aes(x=gene_name, y=cell, fill=geneOCR_table_type)) +
            geom_tile(color="grey50", size=1) +
            scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#42f5ce","#ef42f5"),na.value ="white") +
            theme(
                panel.background = element_blank(),
                axis.title.x=element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_text(angle = 90, color = "grey50", size=7),
                axis.title.y=element_blank()
            ) +
            labs(fill="Gene predicted By")

        if (gene_n >= 20){
            p = p + ggtitle("Most cell variable Gene Top 20")
            gene_n=length(unlist(df %>% distinct(gene_name)))
        }
    }

    list(p=p,width=gene_n,height=cell_n)

}

plot_ocr_gene <- function(tbl_out, plot_type){
    df = tbl_out %>%
        mutate(called=ifelse(called >= 0.02, 1, 0)) %>%
        group_by(chr, start, end, id, cell, gene_id, gene_name) %>%
        summarise(
            called = sum(called),
            geneOCR_table_type=paste(unique(geneOCR_table_type), collapse = ",")
        ) %>%
        ungroup()

    if (plot_type=="heatmap"){
        p = df %>%
            filter(is.na(gene_id)) %>%
            arrange(chr, start) %>%
            mutate(ocr=glue::glue("{chr}:{start}-{end}")) %>%
            mutate(ocr=factor(ocr, levels=unique(ocr))) %>%
            mutate(geneOCR_table_type=ifelse(called==0,NA,geneOCR_table_type)) %>%
            ggplot(aes(x=ocr, y=cell, fill=geneOCR_table_type)) +
            geom_tile(color="grey50", size=1) +
            scale_fill_manual(values=c("#d53e4f","#f46d43","#fdae61","#fee08b","#e6f598","#42f5ce","#ef42f5"),na.value ="white") +
            theme(
                panel.background = element_blank(),
                axis.title.x=element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_text(angle = 90, color = "grey50", size=7),
                axis.title.y=element_blank()
            ) +
            labs(fill="OCR predicted By")
    } else if(plot_type=="arc"){
        Intgraph = graph_from_data_frame(
            df %>%
                filter(is.na(gene_id)) %>%
                mutate(
                    from=df %>%
                        distinct(gene_name) %>%
                        filter(!is.na(gene_name)) %>%
                        pull(gene_name)
                ) %>%
                arrange(chr, start) %>%
                mutate(to=glue::glue("{chr}:{start}-{end}")) %>%
                filter(called >= 0.02) %>%
                dplyr::select(from, to, cell, geneOCR_table_type),
            vertices=df %>%
                mutate(name=glue::glue("{chr}:{start}-{end}")) %>%
                mutate(size=end-start) %>%
                mutate(pos=(start+end)/2) %>%
                mutate(type=ifelse(
                    is.na(id), "promoter", "OCR"
                )) %>%
                dplyr::select(name, size, pos, type) %>%
                distinct() %>%
                bind_rows(
                    df %>%
                        filter(!is.na(gene_id)) %>%
                        dplyr::slice(which.min(end)) %>%
                        distinct(gene_name, end) %>%
                        dplyr::rename(name=gene_name, pos=end) %>%
                        mutate(type="gene", size=0)
                ) %>%
                arrange(pos) %>%
                mutate(name=factor(name, levels=unique(name))),
            directed=F
        )
        p = ggraph(Intgraph, layout="linear") +
            facet_edges(.~cell, drop = F, ncol = 1, strip.position="left") +
            geom_edge_arc(aes(edge_color=geneOCR_table_type), edge_width=0.3, fold=TRUE, show.legend = T) +
            geom_node_point(aes(size=size, color=type), alpha=0.5) +
            geom_node_text(aes(label=name), angle=45, nudge_y = 0, size=1) +
            # scale_x_discrete(labels=as.character(unlist(ocrs %>% select(-fpkm, -cell) %>% distinct() %>% select(dist)))) +
            theme(
                axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
            ) +
            labs(edge_color="Interaction Detected by", size="Genomic Region (bp)", color = "Annotation")
    }
    cell_n = df %>%
        filter(!is.na(cell))%>%
        distinct(cell) %>%
        dplyr::count() %>%
        pull(n)
    list(p=p, width=10, height=cell_n)
}
############################## ui ##############################
ui <- dashboardPage(
    ##### header ####
    dashboardHeader(
        title = "Gene Cis-element Interactome in Pancreas Hi-C (GCIPH) v0.1",
        titleWidth = 450
    ),
    ##### sidebar #####
    dashboardSidebar(
        ##### select cells and Hi-C types #####
        # pickerInput (input$cellInput)
        pickerInput("cellInput","Cells",
                    choices=c("Acinar","Alpha","Beta"),
                    options = list(`actions-box` = TRUE), multiple = T
        ),

        # pickerInput (input$hicInput)
        pickerInput("hicInput","Hi-C CRE-gene prediction models",
                    choices=c("Binary","ABC model","Promoter"),
                    options = list(`actions-box` = TRUE), multiple = T
        ),

        # radioButtons (input$searchItem)
        radioButtons(
            "searchItem", "Please select search item",
            choices = c("SNP", "GENE"), selected = NULL
        ),

        # menu_ui
        uiOutput("menu_ui"),

        ##### action #####
        # button (input$button)
        actionButton(
            "button",label = "search"
        )

    ),
    #### body ####
    dashboardBody(
        shinyjs::useShinyjs(),
        fluidRow(
            # console (output$text in shinyjs::html("text", ""))
            box(
                title = "console", width=12, solidHeader = T, collapsible = TRUE,
                textOutput("text")
            ),

            # plot (output$plot)
            box(
                title = "plot", width=12, solidHeader = TRUE,
                plotOutput("plot"),
                downloadButton("downloadPlot","Download plot")
            ),

            # plot (output$table)
            box(
                title = "table", width=12,
                dataTableOutput("table"),
                downloadButton('downloadTable', 'Download table')
            )
        )
    )
)

############################## server ##############################
server <- function(input, output, session) {

    ##### uiOutput #####
    # menu_ui
    output$menu_ui <- renderUI({
        if(input$searchItem == "SNP"){
            tagList(
                # snp input format widget (input$precal)
                radioButtons(
                    "precal",
                    label = "Select input format",
                    choices = c("Single SNP", "SNiPA output file"),
                    selected = "SNiPA output file"
                ),

                # snp input ui (output$snp_ui)
                # either show widgets for LDlink, or file upload for SNiPA
                uiOutput("snp_ui"),

                # gene filter checkboxInput (input$gene_filter)
                checkboxInput(
                    "gene_filter", label = "Protein coding only",
                    value = T
                ),

                # gene open filter checkboxInput (input$open_gene_only)
                checkboxInput(
                    "open_gene_only", label = "Open gene only",
                    value = T
                )
            )
        }else if(input$searchItem == "GENE"){
            updateSelectizeInput(session, 'gene_name', choices = genenames, server = TRUE)

            tagList(
                # gene_name search selectizeInput (input$gene_name)
                selectizeInput(
                    inputId = "gene_name", label = "Please select a gene name:",
                    choices = NULL, # using updateSelectizeInput for quick load long list of choice
                    multiple = F
                ),

                # uiOutput based on gene_name (output$gene_id_ui)
                uiOutput("gene_id_ui"),

                # gene-OCR plot type (input$plot_type)
                radioButtons(
                    "plot_type", label = "plot type",
                    choices=c("arc","heatmap"),
                    selected = "arc", inline = T
                )
            )
        }
    })


    # snp_ui
    output$snp_ui <- renderUI({
        if(input$precal == "Single SNP"){
            # multiple widgets required here, so tagList is required
            tagList(
                # textInput (input$sentinel)
                textInput(
                    'sentinel', label = "Rs number or snp position",
                    value = "", width = NULL,placeholder = NULL
                ),

                p("   eg:rs10797432,chr1:2501338"),
                # checkboxInput (input$proxy_search)
                checkboxInput(
                    "proxy_search", label = "Search proxy",
                    value = T
                ),
                # ui (output$proxy_ui)
                uiOutput("proxy_ui")
            )

        }else if(input$precal == "SNiPA output file"){
            # fileInput (input$proxy_file)
            fileInput(
                "proxy_file", label = "Upload SNiPA output",
                multiple = F
            )
        }
    })

    # proxy_ui
    output$proxy_ui <- renderUI({
        if(input$proxy_search==T){
            tagList(
                h4(   a(href = 'https://ldlink.nci.nih.gov/?tab=apiaccess', 'apply LDlink API token')),
                # textInput (input$token)
                textInput(
                    'token', label = "LDlink API token",
                    value = "", width = NULL,placeholder = NULL
                ),
                # radioButtons (input$pop)
                radioButtons(
                    "pop", label = "Population",
                    choices = c("AFR", "AMR", "EAS", "EUR", "SAS"), selected = "EUR",inline=T
                ),
                # sliderInput (input$r2)
                sliderInput(
                    "r2", label = "Minimun Rsquare",
                    min=0.1, max=1, value=0.8, step=0.05
                ),
                # sliderInput (input$window_size)
                sliderInput(
                    "window_size", label = "Window size",
                    min=1e5, max=1e6, value=5e5, step=1e5
                )
            )
        }
    })

    # gene_id_ui
    output$gene_id_ui <- renderUI({
        ## selectInput (input$gene_id)
        selectInput(
            "gene_id", label = "Please select a gene id:",
            choices = gene_type %>% filter(gene_name==input$gene_name) %>% pull(gene_id)
        )
    })

    ##### generate table #####
    activeData <- eventReactive(
        input$button,
        {
            withCallingHandlers({
                shinyjs::html("text", "")

                validate(
                    need(
                        !is.null(input$cellInput),
                        "Please select at least one cell"
                    )
                )

                validate(
                    need(
                        !is.null(input$hicInput),
                        "Please select at least one prediction model"
                    )
                )

                if(input$searchItem == "SNP"){
                    ###### SNP tables #######
                    ###### SNP proxy #######
                    if(input$precal=="SNiPA output file"){
                        # validate snipa file is there
                        validate(
                            need(
                                !is.null(input$proxy_file$datapath),
                                "Please load SNiPA output file"
                            )
                        )
                        proxies = read_snipa(input$proxy_file$datapath)
                    }else if(input$precal=="Single SNP"){
                        # validate sentinel format
                        validate(
                            need(
                                (grepl("rs",input$sentinel) | grepl(":",input$sentinel)),
                                 "invalid snp input. Please input rs number or pos (chr:pos) format"
                            )
                        )

                        if (input$proxy_search){
                            # validate LDlink token
                            validate(
                                need(
                                    input$token!="",
                                    "require LDlink API. request API token at https://ldlink.nci.nih.gov/?tab=apiaccess"
                                )
                            )
                        }

                        proxies = cal_proxy(
                            input$proxy_search,
                            input$sentinel,
                            input$window_size,
                            input$pop,
                            input$r2,
                            input$token
                        )
                    }

                    ###### v2g table #######
                    # validate whether proxies table has error
                    validate(
                        need(
                            is.null(proxies$error),
                            proxies$error
                        )
                    )
                    # v2g
                    message("Perform variant to gene mapping...<br>")
                    tbl_out =  do.call(
                        "bind_rows",
                        lapply(
                            input$hicInput, # geneOCR_table_type
                            v2g_mapping_batch,
                            proxy_df = proxies$proxies,
                            cells = input$cellInput
                        )
                    )

                    message("<br>Add query-proxy information...<br>")
                    tbl_out = tbl_out %>%
                            left_join(
                                proxies$proxies %>%
                                dplyr::select(query, proxy, R.squared)
                    )


                    # validate initial v2g is not empty
                    validate(
                        need(
                            nrow(tbl_out)!=0,
                            "no genes were found for given snp(s)"
                        )
                    )

                    # apply gene filter
                    if (input$gene_filter==T){
                        message("<br>Filter genes with protein coding only...<br>")
                        tbl_out = tbl_out %>%
                                semi_join(
                                    gene_type %>%
                                        filter(gene_type=="protein_coding")
                                )

                    }

                    if (input$open_gene_only==T){
                        message("<br>Filter genes with open promoter only...<br>")
                            tbl_out = tbl_out %>%
                                semi_join(
                                    refOCR_cell_simple %>%
                                        filter(!is.na(gene_ids)) %>%
                                        separate_rows(gene_ids, sep="\\|") %>%
                                        dplyr::rename(gene_id=gene_ids)
                                )

                    }

                    # validate filtered v2g is not empty
                    validate(
                        need(
                            nrow(tbl_out)!=0,
                            "no genes were found for given snp(s) after gene filter"
                        )
                    )

                    # generate table for table render
                    message("<br>Generate table ...<br>")
                    tbl_render = tbl_out %>%
                        mutate_at(c("start","end","ocr_start","ocr_end"), function(x){gsub(" ","", format(x,scientific=F))}) %>%
                        mutate(proxy_pos=glue::glue("{chr}:{end}")) %>%
                        left_join(
                            refOCR_cell_simple %>%
                                dplyr::select(id, cell, ocr_enhancer_activity = activity_base, ocr_accessibility_RPKM=ATAC.RPKM)
                        ) %>%
                        mutate(ocr=glue::glue("{ocr_chr}:{ocr_start}-{ocr_end}")) %>%
                        dplyr::select(
                            query, proxy, proxy_pos, R2=R.squared,
                            ocr, ocr_accessibility_RPKM, ocr_enhancer_activity,
                            gene_name, gene_id,
                            cell, called, gene_detected_by = geneOCR_table_type, hicLoops
                        ) %>%
                        filter(called >= 0.02)

                    # generate plot for plot render
                    message("<br>Generate plot ...<br>")
                    p_render = plot_snp_gene(tbl_out)

                }else if(input$searchItem == "GENE"){
                    ##### gene tables ######
                    # validate gene_id
                    validate(
                        need(
                            !is.null(input$gene_id),
                            "Please select corresponding gene_id"
                        )
                    )

                    # generate table
                    message("Search for connected OCRs... <br>")
                    tbl_out =  do.call(
                        "bind_rows",
                        lapply(
                            input$hicInput, # geneOCR_table_type
                            search_geneOCR,
                            gene_id=input$gene_id,
                            cells = input$cellInput
                        )
                    )

                    # validate OCR found for given id
                    validate(
                        need(
                            nrow(tbl_out)!=0,
                            "no OCRs connection were found for given genes"
                        )
                    )

                    # generate table for table render
                    message("<br>Generate table ...<br>")
                    tbl_render = tbl_out %>%
                        mutate_at(c("start","end"), function(x){gsub(" ","", format(x,scientific=F))}) %>%
                        left_join(
                            refOCR_cell_simple %>%
                                dplyr::select(id, cell, ocr_enhancer_activity = activity_base, ocr_accessibility_RPKM=ATAC.RPKM)
                        ) %>%
                        mutate(ocr=glue::glue("{chr}:{start}-{end}")) %>%
                        mutate(query_gene=input$gene_name, query_gene_id=input$gene_id) %>%
                        dplyr::select(
                            query_gene, query_gene_id,
                            ocr, ocr_accessibility_RPKM, ocr_enhancer_activity,
                            cell, called, ocr_detected_by = geneOCR_table_type, hicLoops
                        ) %>%
                        filter(called >= 0.02)

                    # generate plot for plot render
                    message("<br>Generate plot ...<br>")
                    p_render = plot_ocr_gene(tbl_out, plot_type=input$plot_type)
                }

            },
                message = function(m) {
                    shinyjs::html(id = "text", html = m$message, add = TRUE)
                }
            )

            list(tbl_render=tbl_render, p_render=p_render)
        }
    )

    ####### render table output ###########
    output$table <- renderDataTable(
        {
            activeData()[["tbl_render"]]

        },
        escape = FALSE,
        options=list(pageLength = 5)
    )

    output$downloadTable <- downloadHandler(
        filename = function() {
            if (input$searchItem == "SNP") {
                "snpSearch_results.csv"
            }else if (input$searchItem == "GENE") {
                "geneSearch_results.csv"
            }
        },

        content = function(file){
            write.csv(
                activeData()[["tbl_render"]],
                file=file,
                row.names=F
            )
        }
    )
    ######## generate plot output ###########
    output$plot <- renderPlot({
        activeData()[["p_render"]][["p"]]
    }
    )

    output$downloadPlot <- downloadHandler(
        filename = function() {
            if (input$searchItem == "SNP") {
                "snpSearch_results.pdf"
            }else if (input$searchItem == "GENE") {
                "geneSearch_results.pdf"
            }
        },

        content = function(file) {
            pResult = activeData()[["p_render"]]
            ggsave(file,pResult[["p"]],width=pResult[["width"]],height=pResult[["height"]])
        }
    )
}

############################## run ##############################

shinyApp(ui = ui, server = server)
