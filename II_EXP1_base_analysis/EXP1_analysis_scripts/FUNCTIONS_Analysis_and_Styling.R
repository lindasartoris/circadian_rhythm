########################################################################
##########  ADRIANO EXP1: ANALYSIS AND STYLING FUNCTIONS ###############



########## STATS FUNCTIONS ###############
require(report)
require(sjPlot)
require(lme4)
require(car)

# function to test normality of residuals
test_norm <- function(resids) {
  print("Testing normality")
  if (length(resids) <= 300) {
    print("Fewer than 300 data points so performing Shapiro-Wilk's test")
    print(shapiro.test(resids))
    print("below 0.05, the data significantly deviate from a normal distribution")
  } else {
    print("More than 300 data points so using the skewness and kurtosis
approach")
    print("Skewness should be between -3 and +3 (best around zero")
    print(skewness(resids))
    print("")
    print("Excess kurtosis (i.e. absolute kurtosis -3) should be less than 4; ideally around zero")
    print(kurtosis(resids))
  }
}

# function to report a model output
output_lmer <- function(model,show_report = F) {
  print(paste("############### MODEL :",deparse(substitute(m1)),"###############",sep = " "))
  print("------------RESIDUALS NORMALITY------------")
  test_norm(residuals(model))
  print("------------SUMMARY------------")
  print(summary(model))
  print("------------ANOVA------------")
  print(Anova(model))
  print("------------RSQUARED------------")
  print(r.squaredGLMM(model))
  print("------------REPORT------------")
  if (show_report) {
    print(report(model))
  }else{ cat("Set show_report to TRUE to obtain a human readable model report")}
  #tab_model(model)
}

# function to simplify a model if the interaction is not significant
simplify_model <- function(model) {
  # Extract the model formula
  model_formula <- formula(model)
  # Check the significance of the interaction using anova()
  anova_m1 <- as.data.frame(Anova(model))
  print(anova_m1)
  # Find the interaction term in the model formula
  # define a regular expression pattern to match the desired substring
  pattern <- "\\b\\w+\\s*\\*\\s*\\w+\\b"
  # use the sub() function to extract the first match of the pattern
  interaction_term <- sub(paste0(".*(", pattern, ").*"), "\\1", as.character(model_formula)[3])
  interaction_vars <- unlist(strsplit(interaction_term, " * "))
  anova_term <- gsub("\\s*\\*\\s*", ":", interaction_term)
  # If the Anova of the interaction is not significant, simplify the model by removing the interaction
  if (anova_m1[which(rownames(anova_m1)== anova_term),"Pr(>Chisq)"] > 0.05) {
    cat("\n#\nModel interaction NOT significant, simplify\n#\n")
    model_no_interaction_formula <-  as.formula(gsub("\\*", "+", deparse(model_formula)))
    model_no_interaction <- update(model, formula = model_no_interaction_formula)
    #print(summary(m1_no_interaction))
    print(Anova(model_no_interaction))
    return(model_no_interaction)
  } else {
    cat("\n#\nModel interaction significant, don't simplify\n#\n")
    return(model)
  }
}


#check if model has interaction
has_interaction <- function(model) {
  formula_str <- as.character(formula(model))
  return(any(grepl(":", formula_str) | grepl("\\*", formula_str)))
}


# function to perform posthocs
posthoc_list <- list()
interactions_to_explore <- list()
compute_posthocs <- function(model) {
  #warning("this function has only been tested with lmer()")
  warning("How to use: \nA.Run on models without interactions. If interaction is present, run 'simplify_model' first. 
          \nB.If has_interaction= T, paste your variables naming the new var in this format 'VAR1_VAR2'
          \nC. assign the output of this function to 'posthoc_list'
          \nD.Optional: provide an ID_model [i.e. formatted as paste(GROUP,VAR,sep='-')] to later recall the posthoc.\n\n")
  print(paste("model predictor:", paste0(row.names(Anova(model))), sep = " "))
  # check that there are significant outputs
  if (length(row.names(Anova(model)[Anova(model)$"Pr(>Chisq)" < 0.05, ])) == 0) {
    print("there are no significant vars.")
  } else {
    for (SIG.VAR in row.names(Anova(model)[Anova(model)$"Pr(>Chisq)" < 0.05, ])) {
      if (grepl(":", SIG.VAR)) {
        warning(paste0(SIG.VAR, "is an interaction, currently this situation is not handled by the function. Re-run model using pasted variables"))
        #interactions_to_explore <- c(interactions_to_explore, list(paste(GENE,GROUP,SIG.VAR,deparse(substitute(model)), sep = "-") ))
        } else {
        # check if the variable is not numeric . to do so, we need to access the dataframe from the model
        if (!is.numeric(get(gsub("\\[.*", "", as.character(model@call)[3]))[, SIG.VAR])) {
          print(paste0("Performing posthocs for the significant var: ", SIG.VAR))
          arg <- list("Tukey")
          names(arg) <- SIG.VAR
          # glht (mcp) : General linear hypotheses Testing (glht) and multiple comparisons for parametric models (mcp)
          cmp <- do.call(mcp, arg)
          posthoc_SIG.VAR <- summary(glht(model, linfct = cmp), test = adjusted("BH"))
          # Set up a compact letter display of all pair-wise comparisons
          model_means_cld <- cld(posthoc_SIG.VAR)
          # create dataframe usable with ggplot geom_text
          model_means_cld <- as.data.frame(t(t(model_means_cld$mcletters$Letters)))
          # add column name
          model_means_cld$newcol <- NA
          colnames(model_means_cld)[which(names(model_means_cld) == "newcol")] <- SIG.VAR
          colnames(model_means_cld)[which(names(model_means_cld) == "V1")] <- "letters"
          model_means_cld[, SIG.VAR] <- row.names(model_means_cld)
          rownames(model_means_cld) <- NULL
          # if interaction term, split columns
          if (grepl("_", SIG.VAR)) {
            # Split the column into two new columns using "_"
            #SIG.VAR.STRIP <- names(model_means_cld[grep( "_",model_means_cld)])
            model_means_cld[, strsplit(SIG.VAR, "_")[[1]]] <- t(sapply(model_means_cld[,SIG.VAR], function(x) strsplit(x, "_")[[1]]))
          }
          # add to list
          posthoc_list <- c(posthoc_list, list(model_means_cld))
          if (exists("ID_model")) {
            names(posthoc_list)[length(posthoc_list)] <- paste(ID_model,SIG.VAR,  deparse(substitute(model)),sep = "-")
          } else {
            warning("if you provide an ID_model [i.e. formatted as paste(GROUP,VAR,sep='-')] before running the function, it will be easier to later call the posthoc output for plotting")
            names(posthoc_list)[length(posthoc_list)] <- paste(SIG.VAR,deparse(substitute(model)), sep = "-")
          }
          print(paste(deparse(substitute(model)), SIG.VAR, sep = "_"))
          print(model_means_cld)
        } # SIG.VAR LOOP
      } # check if is an interaction
    } # check if numeric
    warning("call 'posthoc_list' to get posthocs")
  } # if significant vars exist
  return(posthoc_list)
}


# # The calculate_weights function takes as input a categorical variable group and an optional argument ratio, which specifies the desired ratio between the weights of the largest and smallest groups
# calculate_weights <- function(group) {
#   print("function to calculate_weights for imbalanced datasets")
#   #make sure that group is a factor
#   group <- as.factor(group)
#   # Calculate the proportions of each group
#   proportions <- table(group) / length(group)
#   print(proportions)
#   # Calculate the weights for each group
#   weights <- numeric(length(group))
#   
#   for (g in levels(group)) {
#     weights[group == g] <- 1 / proportions[g]
#   }
#   # Normalize the weights so that they sum up to 1
#   weights <- weights / sum(unique(weights))
#   return(weights)
# }

## pretty print the model selection output
nice_print_model_sel <- function(model_output) {
  # clean output
  sel.table <- round(as.data.frame(model_output)[-c(1:5)], 3)
  # number of parameters (df) should be K
  names(sel.table)[1] <- "K"
  sel.table$Model <- rownames(sel.table)
  rownames(sel.table) <- NULL
  # replace Model name with formulas
  for (i in 1:nrow(sel.table)) sel.table$formula[i] <- as.character(formula(get(sel.table$Model[i])))[3]
  return(sel.table)
}

# convert significance levels to stars
add_star <- function(p) {
  if (p<0.001) {
    return('***')
  } else if (p<0.01) {
    return('**')
  } else if (p<0.05) {
    return('*')
  } else {
    return('ns')
  }
}

#Box-cox transformation
Box_Cox <- function(x) {
  require(MASS)
  bc <- boxcox(x ~ 1, plotit = FALSE)
  #computes the log-likelihood for a range of lambda values and returns the lambda that maximizes the log-likelihood
  lambda <- bc$x[which.max(bc$y)]
  (x^lambda - 1) / lambda
}


# Define transformations
transformations <- list(
  "Original" = function(x) x,
  "Log" = function(x) log(x - min(x) + 1),
  "Square_Root" = function(x) sqrt(x - min(x)),
  "Box_Cox" = function(x) {
    require(MASS)
    bc <- boxcox(x ~ 1, plotit = FALSE)
    lambda <- bc$x[which.max(bc$y)]
    (x^lambda - 1) / lambda
  }
)

# # Loop through transformations
# for (trans_name in names(transformations)) {
#   # Apply transformation
#   trans_func <- transformations[[trans_name]]
#   transformed_data <- trans_func(sample_data)
#   
# 
#   # Test normality using Shapiro-Wilk test
#   shapiro_test <- shapiro.test(transformed_data)
#   cat(trans_name, "Transformation:\n")
#   cat("Shapiro-Wilk Test p-value:", shapiro_test$p.value, "\n\n")
# }

####################################################
####################################################
####################################################

###########    PLOT SAVING    ###############
require(ggplot2)
SavePrint_plot <- function(plot_obj, plot_name, dataset_name, save_dir, plot_size = c(7, 4), dpi = 300, font_size = 30) {
  # Create the directory if it doesn't exist
  if (!dir.exists(save_dir_plots)) {
    dir.create(save_dir_plots, recursive = TRUE)
  }
  # Modify the plot object to adjust the font size for jpg
  plot_obj_jpg <- plot_obj + theme(text = element_text(size = font_size, lineheight = .3))
  # Check if the directory is writable
  if (!file.access(save_dir_plots, 2)) {
    # Save plot as png
    ggsave(paste0(save_dir_plots, dataset_name, "_", plot_name, "_", Sys.Date(), ".png"), plot = plot_obj_jpg, width = plot_size[1], height = plot_size[2], dpi = dpi)
    # Save plot as pdf
    ggsave(paste0(save_dir_plots, dataset_name, "_", plot_name, "_", Sys.Date(), ".pdf"), plot = plot_obj, width = plot_size[1], height = plot_size[2])
    # Print the plot to the currently open device (the cumulative PDF file)
    print(plot_obj)
  } else {
    cat("Error: The directory is not writable.")
  }
}

########## STYLING FUNCTIONS ###############
require(RColorBrewer)
require(shades)
require(colorspace)
require(plotwidgets)
require(ggplot2)
require(ggnewscale)
require(extrafont)
#library(Cairo) # to ensure working export to pdf with non-standard fonts
library(showtext)
font_add_google("Crimson Text", "crimson", db_cache = FALSE)
showtext_auto()#must be called to indicate that showtext is going to be automatically invoked to draw text whenever a plot is created.


# Import system fonts
#font_import()
# font_import(paths = "/usr/share/texmf/fonts/opentype/public/tex-gyre", pattern = "texgyretermes")
# # Set the default font family
# loadfonts()

#font_import(pattern = "Liberation", prompt= FALSE)
#loadfonts(device = "pdf", quiet = TRUE)

# #Create a custom color scale FOR COLONIES + treatments
# FullPal <- scales::viridis_pal(option = "D")(20)
# myColorsSmall  <- tail(FullPal,5)
# myColorsLarge  <- head(FullPal,5)
# Cols_treatments <- c("#440154FF","#FDE725FF") #"#31688EFF"
# myColors      <- c(myColorsLarge,myColorsSmall, Cols_treatments)
# names(myColors) <- c("R3BP","R5BP","R7BP","R8BP","R12BP","R1SP", "R2SP", "R5SP", "R7SP","R11SP","Big Pathogen","Small Pathogen")
# colScale <- scale_colour_manual(name = "Colony",values = myColors,drop=TRUE)
# fillScale <- scale_fill_manual(name = "Colony",values = myColors,drop=TRUE)

#### CREATE CONSISTENT COLORING FOR ALL THE PLOTTING
# GENERATE 4 MAIN COLORS FOR THE 4 TREATMENTS BS,BP,SS,SP + SHADES FOR THE RESPECTIVE COLONIES

# Create a color palette with 4 colors as distant as possible
colors_full <- scales::viridis_pal(option = "D")(100)
# Create a list to store the subsets of colors
color_subsets <- list()
Shades <- list()
# Loop over the 4 colors to get shades of the colour in a +5, -5 range
for (i in c(10, 30, 70, 90)) { # BE CAREFUL: IF CHANGING THIS, ALSO CHANGE Treat_colors
  color_ramp <- colorRampPalette(colors_full[(i-5):(i+5)])
  #color_ramp <- colorRampPalette(colors_full[(i-1):(i+1)])
  color_subsets[[i]] <- color_ramp(12)
}


Treat_colors <- structure(list(Shades = c(colors_full[10],colors_full[30],colors_full[70],colors_full[90]), 
                               Treatment = c("Big Pathogen", "Big Sham", "Small Pathogen", "Small Sham")),
                          row.names = c(NA, 4L), class = "data.frame")

#show_col(c(colors_full[10],colors_full[30],colors_full[70],colors_full[90]))


#clean the output
color_subsets <- color_subsets[lapply(color_subsets,length)>0]

#Darken the colours progressively to increase contrast among the colonies of the same treatment
for (i in 1:4) {
# Define the color gradient
color_shades <- color_subsets[[i]]

# Convert the colors from hexadecimal to HSL (hue, saturation, lightness)

colors_lightGrad <- c()
# Decrease the lightness of each color by an increasing amount
lightness_decrease <- rev(seq(from = 0, to = 0.2, length.out = length(color_shades)))
lightness_increase <- seq(from = 0, to = 0.2, length.out = length(color_shades))

for (j in 1:length(color_shades)) {
  hsl_colors <- col2hsl(color_shades[j])
  hsl_colors[3] <- hsl_colors[3] - lightness_decrease[j]
  hsl_colors[3] <- hsl_colors[3] + lightness_increase[j]
  colors_lightGrad <- c(colors_lightGrad,hsl2col(hsl_colors))
}

Shades[[i]] <- colors_lightGrad
}

# #inspect output
# par(mfrow = c(2, 2))
# for (i in 1:4) {
#   plot(1:12, 1:12,
#        col = Shades[[i]], # color_subsets
#        pch = 19,
#        cex = 5,
#        xaxt = "n",
#        yaxt = "n",
#        xlab = "",
#        ylab = "") 
# }

### ADD THE METADATA
meta.data <- read.table(paste("/home/cf19810/Documents/scriptsR/EXP1_base_analysis/EXP_summary_data/Metadata_Exp1_2021_2023-02-27.txt", sep = ""), header = T, stringsAsFactors = F, sep = ",")

# create groups to assign the colours
Cols <- list()
# divide each size_treat into a list element with its colonies inside
for (i in 1:length(unique(meta.data$size_treat))) {
  treatment <- unique(meta.data$size_treat)[i]
  treatment_vector <- unique(meta.data$REP_treat[grepl(treatment, meta.data$REP_treat)])
  Cols[[i]] <- treatment_vector
  names(Cols)[i] <- treatment
}
#name list elements according to the favoured pattern (the colour order I have been using since the first plots)
names(Shades)[1] <- "BP"
names(Shades)[2] <- "BS"
names(Shades)[3] <- "SP"
names(Shades)[4] <- "SS"

# dput((list_of_vectors))
# dput((Shades))

# Create an empty dataframe to store the results
colour_palette <- data.frame()

# bind together colours, REP_treats and treatments
# Loop over the list "Cols" to create the dataframe
for (group in names(Cols)) {
  group_cols <- Cols[[group]]
  group_shades <- Shades[[group]]
  # Take a random subset of the shades that matches the length of the cols
  rand_shades <- sample(group_shades, length(group_cols))
  # Create a dataframe with two columns: "Cols" and "Shades"
  group_colour_palette <- data.frame(Cols = group_cols, Shades = rand_shades, Treatment = group)
  # Append the current group dataframe to the overall dataframe
  colour_palette <- rbind(colour_palette, group_colour_palette)
}

# #visualise output
# ggplot(colour_palette, aes(x=Treatment, y=Shades, fill=Shades)) + 
#   geom_tile(colour="white", size=0.5) + 
#   scale_fill_identity() + 
#   theme_void() + 
#   labs(title="Colours by Treatment", x="Treatment", y="Shades") + 
#   theme(axis.text.y=element_text(angle=0, hjust=1)) +
#   facet_wrap(~Treatment, scales = "free_y")

myColors_Colony <- colour_palette$Shades
names(myColors_Colony) <- colour_palette$Cols
myColors_Treatment <- Treat_colors$Shades
names(myColors_Treatment) <- Treat_colors$Treatment

#### DEFINE THE FILL AND COLOR AS STANDALONE

#COLOR
# geom_point(aes(color = Sample)) +
colScale_Colony <- scale_colour_manual(name = "Colony", values = myColors_Colony, drop = TRUE)
#new_scale_color() +
#new_scale_fill() +
# geom_line(aes(color = Sample)) +
colScale_Treatment <- scale_color_manual(name = "Treatment", values = myColors_Treatment) #for lines
colScale_TREATMENT <- scale_color_manual(name = "TREATMENT", values = myColors_Treatment) #for lines
####FILL
# geom_point(aes(color = Sample)) +
colFill_Colony <- scale_fill_manual(name = "Colony", values = myColors_Colony, drop = TRUE)
#new_scale_color() +
#new_scale_fill() +
# geom_line(aes(color = Sample)) +
colFill_Treatment <- scale_fill_manual(name = "Treatment", values = myColors_Treatment) #for lines
colFill_TREATMENT <- scale_fill_manual(name = "TREATMENT", values = myColors_Treatment) #for lines

#### DEFINE REMAINING PLOT STYLE
# ggplot PLOT STYLE
STYLE <- list(
  #colScale, fillScale,
  theme_bw(),
  theme( panel.grid.minor = element_blank(),text=element_text(family="Liberation Serif")),
  scale_x_discrete(labels = function(x) str_wrap(x, width = 4)) # wrap lables when long
)

STYLE_CONT <- list(
  #colScale, fillScale,
  theme_bw(),
  theme(panel.grid.minor = element_blank(),text=element_text(family="Liberation Serif"))
)


remove <- c("color_ramp", "color_shades", "color_subsets", "colors_full", 
            "colors_lightGrad", "colour_palette", "Cols",
            "group", "group_colour_palette", 
            "group_cols", "group_shades", "hsl_colors", "i", "j", "lightness_decrease", 
            "lightness_increase", "meta.data", "myColors_Colony", 
            "rand_shades", 
            "Shades")
# cleaning
rm(list = ls()[which(ls() %in% remove)])
gc()

# #COLOUR SCALES TAKE NAMED VECTORS AS INPUTS
# myColors <- c(Treat_colors$Shades,colour_palette$Shades)
# names(myColors) <- c(Treat_colors$Treatment,colour_palette$Cols)
# colScale <- scale_colour_manual(values = myColors, drop = TRUE) #name = "Colony", 
# fillScale <- scale_fill_manual(values = myColors, drop = TRUE)