repola <- function(D0, h0, n.year, needle_mas = needle_mass, ste = site, params = parameters, reps = repo) {
  if ((ste %in% c("Hyde", "Lettosuo")) == F) {stop("Unknown site: Please pick between Hyde and Lettosuo")}

  reps[c("b0.repo")]<--6.303
  reps[c("b1.repo")]<-14.472
  reps[c("b2.repo")]<--3.976
  reps[c("uk.repo")]<-0.109
  reps[c("eki.repo")]<-0.118

  diameter<-D0*100
  height<-h0

  dski.repo<-2+1.25*diameter

  needle_mas[n.year]<-exp(reps[c("b0.repo")]+reps[c("b1.repo")]*dski.repo/(dski.repo+6)+reps[c("b2.repo")]*height/(height+1))
  needle_mas <- unlist(needle_mas)
  m.N.tot<-needle_mas[n.year]*params[c("carbon_share"),c(ste)]		# kg C / tree
  m.N<-m.N.tot/(params[c("n_age"), c(ste)]*params[c("n_lenght"),c(ste)])				# youngest needles (kg C / mm needle)

  m.R.tot<-m.N.tot*0.5					# kg C / tree

  repola_p <- list(needle_mas, m.N.tot, m.N, m.R.tot)
  names(repola_p) <- c("needle_mass", "m.N.tot", "m.N", "m.R.tot")

  return(repola_p)
}
