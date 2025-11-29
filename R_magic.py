# %%
import rpy2.robjects as robjects
import rpy2.robjects.lib.ggplot2 as ggplot2

from Interfacing_R import seq_data

# %load_ext rpy2.ipython

# %% language="R"
# seq.data <- read.delim('sequence.index', header=TRUE, stringsAsFactors=FALSE)
# seq.data$READ_COUNT <- as.integer(seq.data$READ_COUNT)
# seq.data$BASE_COUNT <- as.integer(seq.data$BASE_COUNT)

# %%
# seq_data = %R seq.data
print(type(seq_data))  #pandas dataframe???

# %%
my_col = list(seq_data.colnames).index("CENTER_NAME")
seq_data['CENTER_NAME'] = seq_data['CENTER_NAME'].apply(lambda x: x.upper())

# %%
# %R -i seq_data
# %R print(colnames(seq_data))

# %% language="R"
# seq_data <- seq_data[seq_data$WITHDRAWN==0, ]
# seq_data$POPULATION <- as.factor(seq_data$POPULATION)

# %% language="R"
# bar <- ggplot(seq_data) +  aes(factor(CENTER_NAME)) + geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# print(bar)

# %% language="R"
# seq_data$POPULATION <- as.factor(seq_data$POPULATION)
# yri_ceu <- seq_data[seq_data$POPULATION %in% c("YRI", "CEU") & seq_data$BASE_COUNT < 2E9 & seq_data$READ_COUNT < 3E7, ]

# %% language="R"
# scatter <- ggplot(yri_ceu, aes(x=BASE_COUNT, y=READ_COUNT, col=factor(ANALYSIS_GROUP), shape=POPULATION)) + geom_point()
# print(scatter)

# %% language="R"
# library(gridExtra)
# library(grid)
# g <- grid.arrange(bar, scatter, ncol=1)
# g

# %% language="R"
# png('fig.png')
# g
# dev.off()