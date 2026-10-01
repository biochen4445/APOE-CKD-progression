# Purpose: eGFR equations (mL/min/1.73 m^2). Creatinine in mg/dL, age in years, sex "M"/"F".

# CKD-EPI 2021 race-free creatinine equation (Inker et al., NEJM 2021)
egfr_ckdepi2021 <- function(scr, age, sex) {
  female <- sex == "F"
  kappa  <- ifelse(female, 0.7, 0.9)
  alpha  <- ifelse(female, -0.241, -0.302)
  142 * pmin(scr / kappa, 1)^alpha * pmax(scr / kappa, 1)^(-1.200) *
    0.9938^age * ifelse(female, 1.012, 1)
}

# 4-variable MDRD without race term. constant = 175 (IDMS-traceable) or 186.3 (JAMA 2005)
egfr_mdrd <- function(scr, age, sex, constant = 175) {
  constant * scr^(-1.154) * age^(-0.203) * ifelse(sex == "F", 0.742, 1)
}

egfr_calc <- function(scr, age, sex, equation = c("ckdepi2021", "mdrd"), mdrd_constant = 175) {
  equation <- match.arg(equation)
  switch(equation,
         ckdepi2021 = egfr_ckdepi2021(scr, age, sex),
         mdrd       = egfr_mdrd(scr, age, sex, mdrd_constant))
}
