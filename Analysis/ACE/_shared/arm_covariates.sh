#!/bin/bash
# Single source of truth for the non-ANCOVA male ACE AD-model arms (shell side).
# Mirrors arm_covariates.py / arm_covariates.R. Sourced by the per-arm launchers.
# All arms males-only, phenotype tot_adverse_exp, integration derived_batch.

# Arms to iterate
NON_ANCOVA_ARMS=(MaleNoADadj MaleNiaReagan MaleBinaryAD MaleContAD MaleAceByAD)

# Signal-bearing cell types used across all downstream analyses.
# DEG/GSEA/TF use label "Inh"; SCENIC/WGCNA h5ad input uses file "broad_Inh".
SIGNAL_CELLTYPES=(In-PV_Basket Ast Mic Oli OPC Inh)

# arm -> space-separated AD covariates added to "age_death pmi"
arm_ad_covars() {
  case "$1" in
    MaleNoADadj)   echo "" ;;
    MaleNiaReagan) echo "niareagansc" ;;
    MaleBinaryAD)  echo "AD_binary" ;;
    MaleContAD)    echo "amylsqrt tangsqrt" ;;
    MaleAceByAD)   echo "AD_binary" ;;
    *) echo "ERROR_unknown_arm" ;;
  esac
}

# arm -> 1 if it uses an AD_binary:phenotype interaction, else 0
arm_interaction() {
  case "$1" in
    MaleAceByAD) echo 1 ;;
    *) echo 0 ;;
  esac
}

# arm -> DEG .rda object name
arm_deg_obj() {
  case "$1" in
    MaleBinaryAD|MaleAceByAD) echo "res_ace" ;;
    *) echo "res" ;;
  esac
}

# arm -> DEG filename suffix (the ACE-main contrast)
arm_deg_suffix() {
  case "$1" in
    MaleBinaryAD) echo "MaleBinaryAD_ACEmain" ;;
    MaleAceByAD)  echo "MaleAceByAD_ACEmain" ;;
    *) echo "$1" ;;
  esac
}

# Map a cell-type label to its h5ad / SCENIC-dir basename (Inh -> broad_Inh)
ct_to_h5ad() {
  case "$1" in
    Inh) echo "broad_Inh" ;;
    *) echo "$1" ;;
  esac
}
