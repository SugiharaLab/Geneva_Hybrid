require(tidyverse)

########################################################################################
# The original dataset contains the following variables:
# pp = primary production
# chla = concentration in Chlorophyll A
# oxydeep = concentration of oxygen in the hypolimnion of Lake Geneva
# PT = concentration of total phosphorus
# PO4 = concentration in orthophosphate (solve reactive phosphorus)
# Phypso = the total mass of phosphorus present in the water column of Lake Geneva
# thermo = depth of the termocline
# oxycol = oxygen concentration in the water column
# hypoarea = area of the lake for which oxygen concentration is below 2 mg.l-1 (hypoxic)
# surftemp = surface temperature
# rhone = discharge of the Rh?ne river (the main tributary of the lake)
########################################################################################

dictionary_short <- list("DO_b50"="oxydeep",
                         "DOeu"="oxycol",
                         "Ahyp"="hypoxarea",
                         "Pt"="pt",
                         "PO4"="po4",
                         "MP"="phypso",
                         "Zt"="thermo",
                         "T0"="surftemp",
                         "W"="wind_speed",
                         "TAt"="T_air",
                         "QRh"="rhone",
                         "yearday"="yearday",
                         "PP"="pp",
                         "Ch"="chla" )     


# dictionary used for Fig 1D labeling:
dictionary_extended <- list("DO"="oxydeep",
                         "TP_surf"="Ptot_epi",
                         "TP_lake"="Ptot_lake",
                         "SRP_surf"="PO4_epi",
                         "SRP_lake"="PO4_lake",
                         "h_mix"="h_mix",
                         "h_mix_simstrat"="h_mix",
                         "T_surf"="T_surf",
                         "W"="wind_speed",
                         "Tatm"="T_air",
                         "Q"="rhone",
                         "year_sine"="year_sine",
                         "yearday"="yearday",
                         "PP"="pp",
                         "Chl"="chla" )     


# dictionary used for Fig 1D labeling:
# dictionary_expressions <- list("DO"="oxydeep",
#                             "TP[{surf}]"="Ptot_epi",
#                             "TP[{lake}]"="Ptot_lake",
#                             "SRP[{surf}]"="PO4_epi",
#                             "SRP[{lake}]"="PO4_lake",
#                             "h[{mix}]"="Z_thermo",
#                             "h[{mix}]sim"="Z_thermo_model",
#                             "T[{surf}]"="T_surf",
#                             "W"="wind_speed",
#                             "T[{atm}]"="T_air",
#                             "Q"="rhone",
#                             "yearday"="yearday",
#                             "PP"="pp",
#                             "Chl"="chla" )  


dictionary_expressions <- list("DO"="oxydeep",
                               "TP[surf]"="Ptot_epi",
                               "TP[lake]"="Ptot_lake",
                               "SRP[surf]"="PO4_epi",
                               "SRP[lake]"="PO4_lake",
                               "h[mix]"="h_mix",
                               "h[mix]sim"="h_mix_model",
                               "T[surf]"="T_surf",
                               "W"="wind_speed",
                               "T[atm]"="T_air",
                               "Q"="rhone",
                               # "sin(\u03BD)"="year_sine", # on windows machines, UTF-8 encoding probably won't work.
                               "sin(u)"="year_sine",
                               "yearday"="yearday",
                               "PP"="pp",
                               "Chl"="chla" )  

dictionary_extended_to_expressions <- names(dictionary_expressions) %>% set_names(nm = names(dictionary_extended))

f_q_exp <- function(breaks) { parse(text=breaks)}

label_embeddings <- function(embeddings,column_names,dictionary=dictionary_short){
    column_short <- str_replace_all(column_names,set_names(x=names(dictionary),nm=dictionary))

    labels <- map_chr(embeddings,function(embedding){
        embedding <- str_extract_all(embedding,"[:digit:]+")[[1]] %>%
            as.numeric()
        paste0("<",paste(column_short[embedding],collapse=","),">")
    })
    
    # return(as.factor(labels))
    return(labels)
}

label_embeddings_parsable <- function(embeddings,column_names,dictionary=dictionary_short){
    column_short <- str_replace_all(column_names,set_names(x=names(dictionary),nm=dictionary))
    
    labels <- map_chr(embeddings,function(embedding){
        embedding <- str_extract_all(embedding,"[:digit:]+")[[1]] %>%
            as.numeric()
        paste0("\"<\"~",paste(column_short[embedding],collapse="~\",\"~"),"~\">\"")
        # paste0("c~",paste(column_short[embedding],collapse="~\",\"~"),"~c")
        # paste(column_short[embedding],collapse="~\",\"~")
        # column_short[embedding][1]
    })
    
    # return(as.factor(labels))
    return(labels)
}

make_label_parsable <- function(labels) {
    map_chr(labels,function(label){
    embedding <- str_extract_all(label,"[^<>,]+") %>% unlist()
    relabel <- paste0("\"<\"~",paste(embedding,collapse="~\",\"~"),"~\">\"")
    return(relabel)})
}


embedding_int_to_chr <- function(embeddings,column_names){
    
    embeddings <- map(embeddings,function(embedding){
        embedding <- str_extract_all(embedding,"[:digit:]+")[[1]] %>%
            as.numeric()
        return(column_names[embedding])
    })
    
    # return(as.factor(labels))
    return(embeddings)
}


df_dict_plotting <- bind_rows(
    data.frame(data_file="oxydeep", plotting_long="Dissolved Oxygen\nBottom 50m (mg/L)", plotting_short="DO_b50m"),
    data.frame(data_file="T_surf", plotting_long="Lake Surface Temperature (ºC)", plotting_short="T_surf"),
    data.frame(data_file="PO4_epi", plotting_long="Soluble Reactive Phosphorous\nSurface (\u03BCg/L)", plotting_short="SRP_surf"),
    data.frame(data_file="PO4_lake", plotting_long="Soluble Reactive Phosphorous\nLake Averaged (\u03BCg/L)", plotting_short="SRP_lake"),
    data.frame(data_file="Ptot_both",plotting_long="Total Phosphorous\n(\u03BCg/L)", plotting_short="Ptot_both"),
    data.frame(data_file="Ptot_epi",plotting_long="Total Phosphorous\nSurface (\u03BCg/L)", plotting_short="TP_surf"),
    data.frame(data_file="Ptot_lake",plotting_long="Total Phosphorous\nLake Averaged (\u03BCg/L)",plotting_short="TP_lake"),
    data.frame(data_file="rhone",plotting_long="Rhone River Discharge\n(m^3/s)",plotting_short="Q"),
    data.frame(data_file="T_air",plotting_long="Air Temperature (ºC)",plotting_short="T_atm"),
    data.frame(data_file="wind_speed",plotting_long="Wind Speed (m/s)",plotting_short="W"),
    data.frame(data_file="h_mix",plotting_long="Height of\nMixed Layer (m)",plotting_short="h_mix"),
    data.frame(data_file="chl",plotting_long="Chlorophyll-A\n(\u03BCg/L)",plotting_short="Chl-A")
)

L_dict_plotting_long <- df_dict_plotting %>% pull(plotting_long) %>% 
    set_names(df_dict_plotting$data_file)
L_dict_plotting_short <- df_dict_plotting %>% pull(plotting_short) %>%
    set_names(df_dict_plotting$data_file)

L_colors <- bind_rows(
    data.frame(data_file="oxydeep",color="blue"),
    data.frame(data_file="chl",color="green"),
    data.frame(data_file="PO4_epi",color="plum3"),
    data.frame(data_file="PO4_lake",color="plum4"),
    data.frame(data_file="Ptot_epi",color="mediumorchid3"),
    data.frame(data_file="Ptot_lake",color="mediumorchid4")
) %>% {set_names(pull(.,color),pull(.,data_file))}

