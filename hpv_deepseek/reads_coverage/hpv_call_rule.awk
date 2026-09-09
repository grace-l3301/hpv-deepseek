# From reads_coverage_with_ihop_v16_04_02_2029.py (v16), 2026-09-08.

BEGIN { FS = "\t"; res = "NEGATIVE"; if (glen <= 0) glen = 7906; if (tlen <= 0) tlen = 170
        # borderline band: the LOWER, more permissive pair. Unset => no BORDERLINE state.
        if (bordreads == "") bordreads = 0
        if (bordcov   == "") bordcov   = 0
        use_border = (bordreads > 0 || bordcov > 0)
        if (use_border) {
            if (bordreads <= 0) bordreads = minreads
            if (bordcov   <= 0) bordcov   = mincov
        } }
{
    contig = $1; reads = $2 + 0; cov = $3 + 0
    if (contig !~ /[Hh][Pp][Vv]/) next
    n = int(reads / 2)                      # paired records -> templates
    exp_cov = (n > 0) ? 100 * (1 - exp(n * log(1 - tlen / glen))) : 0
    ratio = (exp_cov > 0) ? cov / exp_cov : 0
    pass_r = (reads >= minreads)
    pass_c = (cov   >= mincov)
    verdict = (pass_r && pass_c) ? "PASS" : "fail"
    if (pass_r && pass_c) res = "POSITIVE"
    why = (pass_r ? "reads_ok" : "reads_low") "," (pass_c ? "cov_ok" : "cov_low")
    # Clears the permissive pair but not the strict one => the call is threshold-dependent.
    # Never downgrades a POSITIVE: res is only touched while it is still NEGATIVE.
    if (use_border && !(pass_r && pass_c) && reads >= bordreads && cov >= bordcov) {
        verdict = "borderline"
        why = why ",in_band"
        if (res == "NEGATIVE") res = "BORDERLINE"
    }
    flag = (ratio > 0 && ratio < 0.25 && n >= 5) ? "CLUSTERED_amplicon_like" : "."
    printf "%s\t%s\t%s\t%d\t%.4f\t%.2f\t%.2f\t%.3f\t%s\t%s\t%s\n",
           smp, stage, contig, reads, cov, minreads + 0, mincov + 0, ratio, why, verdict, flag >> audit
}
END { print res }
