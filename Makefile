# ----  R pipeline Makefile ----

R := Rscript
MAIN := R/clustering.R
OUT := project/outputs/figures

# Expected outputs
FIGS := $(OUT)/kmeans_cluster.png \
        $(OUT)/shell_clusters_interative.html \
        $(OUT)/spectral_cluster_1.png\
        $(OUT)/spectral_cluster_0.8.png\
        $(OUT)/spectral_cluster_1.2.png

all: $(FIGS)

$(OUT):
	mkdir -p $(OUT)

$(FIGS): $(MAIN) | $(OUT)
	$(R) $(MAIN)

clean:
	rm -f $(OUT)/*

.PHONY: all clean

