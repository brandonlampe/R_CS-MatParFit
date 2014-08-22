



  Z3  <- rep(0, size(ELC)[2])
  DT.I <- DT
  POW <- rep(0, size(ELC)[2])
  ERROR.OUT <- rep(1, size(ELC)[2])
  ERROR.INN <- rep(1, size(ELC)[2])
  # FU.I calculated with vector of zeros (Z3)
  FU.I <- ifelse(Z3 == EFT, 1, ifelse(Z3 < EFT, exp(BIGD * (1 - Z3 / EFT) ^ 2),
                                      exp(-DELTA * (1 - Z3 / EFT) ^ 2)))
  DZ3.I <- (FU.I - 1) * ESS # derivative calculated

  # derivative integrated to yield orignal values (first iteration)
  Z3.I <- ifelse(TIME > 0, 0.5 * DT * (shift(DZ3.I, 1, "right") + DZ3.I), 0)

  Z3.CHECK <- Z3.I # for debugging
  IT.OUT <- 0
  IT.INN <- 0
  IMAX <- 100
  TOL <- 10 ^ 3 * .Machine$double.eps
  IN.TIME <- proc.time()
  ERROR.OUT <- rep(1, size(ELC)[2])

  repeat{
    Z3.OUT <- Z3.I
    FU.I <- ifelse(Z3.OUT == EFT, {1},{ifelse(Z3.OUT < EFT,{
      exp(BIGD * (1 - Z3.OUT / EFT) ^ 2)},{
      exp(-DELTA * (1 - Z3.OUT / EFT) ^ 2)})})
    DZ3.I <- (FU.I - 1) * ESS
    Z3.I <- ifelse(TIME > 0, 0.5 * DT * (shift(DZ3.I, 1, "right") + DZ3.I), 0)
    ERROR.OUT <- abs(Z3.I - Z3.OUT)
    IT.OUT <- IT.OUT + 1
    print(IT.OUT)
    if(all(ERROR.OUT < TOL)) {
      Z3 <- Z3.I
      break}
    if(any(ERROR.OUT > TOL) & IT.OUT >= IMAX){break}}

  ERROR.INN <- ifelse(ERROR.OUT > TOL, ERROR.OUT, 0)
  Z3.OUT    <- ifelse(ERROR.OUT > TOL, -1, Z3.OUT)
  Z3.INN    <- ifelse(Z.OUT != -1, {0}, {

  })


  foreach()
        repeat{
          POW <- POW + 1        # exponent on divisor
          DT.H <- DT.I/(2^POW)  # divosor
          Z3.INN <- Z3.I
          FU.I <- ifelse(Z3.INN == EFT, {1},{ifelse(Z3.INN < EFT,{
            exp(BIGD * (1 - Z3.INN / EFT) ^ 2)},{
            exp(-DELTA * (1 - Z3.INN / EFT) ^ 2)})})
          Z3.I <- ifelse(TIME > 0, 0.5 * DT.H * (shift(DZ3.I, 1, "right") + DZ3.I), 0)
          ERROR.INN <- abs(Z3.I - Z3.INN)
          IT.INN <- IT.INN + 1  # iterator count on inner loop
          print(ERROR.INN)
          if(ERROR.INN <= TOL || IT.INN >= IMAX){break}}}}}


# ---- SPLIT TIME ----
#     if(ERROR.OUT > TOL){
# #       while((rev(ERROR.INN) < TOL) && (IT.INN < IMAX)){
#         POW <- POW + 1
#         IT.INN <- IT.INN + 1
#         DT.H <- DT.I/(2^POW)
#         Z3 <- Z3.I
#         FU.I <- ifelse(Z3 == EFT, {1},{ifelse(Z3 < EFT,{
#           exp(BIGD * (1 - Z3 / EFT) ^ 2)},{
#             exp(-DELTA * (1 - Z3 / EFT) ^ 2)})})
#         Z3.I <- ifelse(TIME > 0, 0.5 * DT.H * (shift(DZ3.I, 1, "right") + DZ3.I), 0)
#         ERROR.INN <- abs(Z3.I - Z3)
#         print(c("INNER WITH", IT.INN, ERROR.INN))}
# }
