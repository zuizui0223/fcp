#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
  library(rotl)
  library(jsonlite)
})

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag, default=NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(paste("Missing value for", flag))
  args[[i+1]]
}
dataset <- get_arg("--dataset")
outdir <- get_arg("--outdir")
replicates <- as.integer(get_arg("--replicates", "100"))
seed <- as.integer(get_arg("--seed", "20260925"))
if (is.null(dataset) || is.null(outdir)) stop("--dataset and --outdir are required")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

d <- read.csv(dataset, stringsAsFactors=FALSE, check.names=FALSE)
required <- c("species","state","white_fraction","n_classifiable")
if (!all(required %in% names(d))) stop("missing required columns")
d <- d[d$state %in% c("nonwhite_dominant","white_dominant"), , drop=FALSE]
if (nrow(d) < 150) stop(sprintf("endpoint panel below predeclared threshold: %d", nrow(d)))

matches <- tnrs_match_names(d$species, context_name="Land plants", do_approximate_matching=TRUE)
write.csv(matches, file.path(outdir,"tnrs_matches.csv"), row.names=FALSE)
if (!all(c("search_string","ott_id") %in% names(matches))) stop("unexpected TNRS columns")
if (!"is_approximate_match" %in% names(matches)) matches$is_approximate_match <- FALSE
matches$search_norm <- tolower(trimws(as.character(matches$search_string)))
d$query_norm <- tolower(trimws(as.character(d$species)))
input_map <- setNames(d$species, d$query_norm)
matches$input_species <- unname(input_map[matches$search_norm])
matches$ott_id <- suppressWarnings(as.integer(matches$ott_id))
usable <- matches[!is.na(matches$ott_id) & !matches$is_approximate_match & !is.na(matches$input_species), , drop=FALSE]
dup_ott <- usable$ott_id[duplicated(usable$ott_id) | duplicated(usable$ott_id, fromLast=TRUE)]
usable <- usable[!(usable$ott_id %in% dup_ott), , drop=FALSE]
if (nrow(usable) < 150) stop(sprintf("too few exact OpenTree matches: %d", nrow(usable)))

# OpenTree TNRS can return an OTT id that the induced-tree endpoint reports as pruned.
# Handle this outcome-blind taxonomy-availability case mechanically: remove only OTT ids
# explicitly named by the API as pruned, record them, and retry without changing any
# biological state threshold or estimability gate.
pruned_ids <- integer(0)
remaining_ott <- usable$ott_id
repeat {
  ans <- tryCatch(
    list(tree=tol_induced_subtree(ott_ids=remaining_ott, label_format="id"), error=NULL),
    error=function(e) list(tree=NULL,error=conditionMessage(e))
  )
  if (!is.null(ans$tree)) {
    tr <- ans$tree
    break
  }
  msg <- ans$error
  hits <- regmatches(msg, gregexpr("ott[0-9]+", msg))[[1]]
  if (length(hits)==0 || identical(hits,-1)) stop(paste("OpenTree induced_subtree failed:",msg))
  ids <- unique(as.integer(sub("^ott","",hits)))
  ids <- ids[ids %in% remaining_ott]
  if (length(ids)==0) stop(paste("OpenTree induced_subtree failed without removable OTT id:",msg))
  pruned_ids <- unique(c(pruned_ids,ids))
  remaining_ott <- setdiff(remaining_ott,ids)
  if (length(remaining_ott)<150) stop(sprintf("too few OpenTree ids after API-pruned exclusions: %d",length(remaining_ott)))
}
if (length(pruned_ids)) {
  write.csv(data.frame(ott_id=pruned_ids),file.path(outdir,"opentree_api_pruned_ids.csv"),row.names=FALSE)
  usable <- usable[usable$ott_id %in% remaining_ott,,drop=FALSE]
}
write.tree(tr, file=file.path(outdir,"opentree_induced_topology_ott.tre"))
extract_ott <- function(x) as.integer(sub("^ott", "", x))
tip_ott <- vapply(tr$tip.label, extract_ott, integer(1))
map <- setNames(usable$input_species, as.character(usable$ott_id))
tip_names <- unname(map[as.character(tip_ott)])
if (any(is.na(tip_names))) stop("failed to map all tips")
tr$tip.label <- gsub(" ", "_", tip_names, fixed=TRUE)

matched_species <- gsub("_", " ", tr$tip.label, fixed=TRUE)
md <- d[d$species %in% matched_species, , drop=FALSE]
md <- md[match(matched_species, md$species), , drop=FALSE]
rownames(md) <- tr$tip.label
state_counts <- table(md$state)
if (length(state_counts) != 2 || any(state_counts < 50)) {
  stop(paste("state count gate failed:", paste(names(state_counts), state_counts, collapse=", ")))
}

fit_one <- function(tree, states) {
  er <- tryCatch(ace(states, tree, type="discrete", model="ER", method="ML"), error=function(e) NULL)
  ard <- tryCatch(ace(states, tree, type="discrete", model="ARD", method="ML"), error=function(e) NULL)
  if (is.null(er) || is.null(ard)) return(NULL)
  idx <- ard$index.matrix
  if (is.null(idx) || nrow(idx) != 2 || ncol(idx) != 2) return(NULL)
  q_nw <- as.numeric(ard$rates[idx[1,2]])
  q_wn <- as.numeric(ard$rates[idx[2,1]])
  if (!is.finite(q_nw) || !is.finite(q_wn) || q_nw <= 0 || q_wn <= 0) return(NULL)
  ll_er <- as.numeric(er$loglik)
  ll_ard <- as.numeric(ard$loglik)
  lrt <- max(0, 2*(ll_ard-ll_er))
  p <- pchisq(lrt, df=1, lower.tail=FALSE)
  data.frame(q_nonwhite_to_white=q_nw, q_white_to_nonwhite=q_wn,
             rate_ratio=q_nw/q_wn, log_rate_ratio=log(q_nw/q_wn),
             loglik_er=ll_er, loglik_ard=ll_ard, lrt=lrt, lrt_p=p)
}

set.seed(seed)
states <- factor(md$state, levels=c("nonwhite_dominant","white_dominant"))
names(states) <- rownames(md)
rows <- list()
for (r in seq_len(replicates)) {
  tree_r <- if (is.binary.tree(tr)) tr else multi2di(tr, random=TRUE)
  tree_r <- compute.brlen(tree_r, method="Grafen", power=1)
  st <- states[tree_r$tip.label]
  z <- fit_one(tree_r, st)
  if (!is.null(z)) {
    z$replicate <- r
    rows[[length(rows)+1L]] <- z
  }
}
if (length(rows) == 0) stop("no successful transition fits")
raw <- do.call(rbind, rows)
raw <- raw[,c("replicate","q_nonwhite_to_white","q_white_to_nonwhite","rate_ratio","log_rate_ratio","loglik_er","loglik_ard","lrt","lrt_p")]
write.csv(raw, file.path(outdir,"transition_rate_replicates.csv"), row.names=FALSE)

completed <- nrow(raw)
ratio_gt1_fraction <- mean(raw$rate_ratio > 1)
lrt_sig_fraction <- mean(raw$lrt_p < 0.05)
median_ratio <- median(raw$rate_ratio)
median_lrt_p <- median(raw$lrt_p)

status <- "DIRECTIONAL_ASYMMETRY_NOT_SUPPORTED_UNDER_THIS_TEST"
if (length(tr$tip.label) >= 150 && all(state_counts >= 50) && completed >= 80 &&
    median_ratio > 1 && ratio_gt1_fraction >= 0.90 &&
    median_lrt_p < 0.05 && lrt_sig_fraction >= 0.80) {
  status <- "NONWHITE_TO_WHITE_ASYMMETRY_SUPPORTED"
}

# descriptive 0.05 / 0.95 sensitivity using the same tree if enough endpoint tips
sens <- list(status="not_estimable")
if (all(c("white_fraction") %in% names(d))) {
  ds <- d[d$white_fraction <= 0.05 | d$white_fraction >= 0.95, , drop=FALSE]
  ds$state_sens <- ifelse(ds$white_fraction <= 0.05, "nonwhite_dominant", "white_dominant")
  sms <- ds[ds$species %in% matched_species, , drop=FALSE]
  sms <- sms[match(matched_species[matched_species %in% sms$species], sms$species), , drop=FALSE]
  counts_s <- table(sms$state_sens)
  if (length(counts_s)==2 && all(counts_s >= 30)) {
    tr_s <- drop.tip(tr, setdiff(tr$tip.label, gsub(" ","_",sms$species,fixed=TRUE)))
    ss <- factor(sms$state_sens[match(gsub("_"," ",tr_s$tip.label,fixed=TRUE), sms$species)],
                 levels=c("nonwhite_dominant","white_dominant"))
    names(ss) <- tr_s$tip.label
    srows <- list()
    set.seed(seed+1)
    for (r in seq_len(replicates)) {
      tree_r <- if (is.binary.tree(tr_s)) tr_s else multi2di(tr_s, random=TRUE)
      tree_r <- compute.brlen(tree_r, method="Grafen", power=1)
      z <- fit_one(tree_r, ss[tree_r$tip.label])
      if (!is.null(z)) srows[[length(srows)+1L]] <- z
    }
    if (length(srows)>0) {
      sr <- do.call(rbind,srows)
      sens <- list(status="complete", endpoint_counts=as.list(counts_s),
                   completed_replicates=nrow(sr),
                   median_rate_ratio=median(sr$rate_ratio),
                   fraction_rate_ratio_gt1=mean(sr$rate_ratio>1),
                   median_lrt_p=median(sr$lrt_p),
                   fraction_lrt_p_lt_0_05=mean(sr$lrt_p<0.05))
    }
  }
}

summary <- list(
  schema="fcp_white_transition_direction_v1",
  status=status,
  requested_endpoint_species=nrow(d),
  exact_nonapproximate_tnrs_matches=nrow(usable),
  induced_tree_species=length(tr$tip.label),
  state_counts=as.list(state_counts),
  replicates_requested=replicates,
  completed_replicates=completed,
  median_q_nonwhite_to_white=median(raw$q_nonwhite_to_white),
  median_q_white_to_nonwhite=median(raw$q_white_to_nonwhite),
  median_rate_ratio=median_ratio,
  rate_ratio_q025=as.numeric(quantile(raw$rate_ratio,0.025)),
  rate_ratio_q975=as.numeric(quantile(raw$rate_ratio,0.975)),
  fraction_rate_ratio_gt1=ratio_gt1_fraction,
  median_lrt_p=median_lrt_p,
  fraction_lrt_p_lt_0_05=lrt_sig_fraction,
  tree="OpenTree induced topology; random polytomy resolution; Grafen power=1 branch lengths",
  sensitivity_005_095=sens,
  hard_nonclaims=c(
    "does not estimate absolute evolutionary rates",
    "does not establish within-population mutation direction",
    "does not identify anthocyanin loss as the molecular mechanism",
    "does not establish heat or pollinator causation",
    "does not establish the ancestral flower colour of angiosperms"
  )
)
write_json(summary, file.path(outdir,"result.json"), auto_unbox=TRUE, pretty=TRUE)
write.csv(data.frame(species=matched_species,state=md$state,white_fraction=md$white_fraction,
                     n_classifiable=md$n_classifiable),
          file.path(outdir,"matched_endpoint_species.csv"), row.names=FALSE)
print(summary)
