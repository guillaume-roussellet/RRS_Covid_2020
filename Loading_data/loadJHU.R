# ===========================================
# JHU covid data:
library(stringr)
# https://covidtracking.com/api
# =====================================================
# LOAD JHU DATA
# =====================================================
covid_daily <- read.csv("https://covidtracking.com/api/v1/states/daily.csv")

year      <- str_sub(covid_daily$date,1,4)
month     <- str_sub(covid_daily$date,5,6)
day       <- str_sub(covid_daily$date,7,8)
date.aux  <- paste(year,"/",month,"/",day,sep="")

covid_daily$date <- as.Date(date.aux,"%Y/%m/%d")

vector.of.states  <- unique(as.character(covid_daily$state))
indic.excluded    <- which(vector.of.states %in% c("AS","GU","MP","VI","PR"))
vector.of.states  <- vector.of.states[-indic.excluded]

nb.states <- length(vector.of.states)

vector.of.dates  <- as.Date(levels(as.factor(covid_daily$date)))

vector.of.variables <- c("positiveCasesViral","recovered","death")

DATA <- data.frame(date=as.Date(vector.of.dates))
for(i in 1:nb.states){
  State <- vector.of.states[i]
  DATA.aux <- subset(covid_daily,state==State)
  
  for(variable in vector.of.variables){
    eval(parse(text = gsub(" ","",paste("DATA.4.one.state <- data.frame(date=DATA.aux$date,",
                                        variable,".",State,"=DATA.aux$",variable,")", sep=""))))
    
    DATA <- merge(DATA,DATA.4.one.state,by="date",all=TRUE)
  }
}
q = which(DATA[,"date"] == "2020-09-30")
DATA = DATA[1:q,]


DEATH <- data.frame(date=as.Date(vector.of.dates))
for(i in 1:nb.states){
  State <- vector.of.states[i]
  DATA.aux <- subset(covid_daily,state==State)
  
  eval(parse(text = gsub(" ","",paste("DATA.4.one.state <- data.frame(date=DATA.aux$date,death.",State,"=DATA.aux$death)", sep=""))))
  
  DEATH <- merge(DEATH,DATA.4.one.state,by="date",all=TRUE)
}
q = which(DEATH[,"date"] == "2020-09-30")
DEATH = DEATH[1:q,]
DEATH[is.na(DEATH)] <- 0
DEATH[2:dim(DEATH)[1],2:(nb.states+1)] <- DEATH[2:dim(DEATH)[1],2:(nb.states+1)] - DEATH[1:(dim(DEATH)[1]-1),2:(nb.states+1)]

#============================================
# INFECTED AND RECOVERED DATA
#============================================
INFEC <- data.frame(date=as.Date(vector.of.dates))
for(i in 1:nb.states){
  State <- vector.of.states[i]
  DATA.aux <- subset(covid_daily,state==State)
  
  eval(parse(text = gsub(" ","",paste("DATA.4.one.state <- data.frame(date=DATA.aux$date,positiveCasesViral.",
                                      State,"=DATA.aux$positiveCasesViral)", sep=""))))
  
  INFEC <- merge(INFEC,DATA.4.one.state,by="date",all=TRUE)
}
q = which(INFEC[,"date"] == "2020-09-30")
INFEC = INFEC[1:q,]

RECOV <- data.frame(date=as.Date(vector.of.dates))
for(i in 1:nb.states){
  State <- vector.of.states[i]
  DATA.aux <- subset(covid_daily,state==State)
  
  eval(parse(text = gsub(" ","",paste("DATA.4.one.state <- data.frame(date=DATA.aux$date,recovered.",State,"=DATA.aux$recovered)", sep=""))))
  
  RECOV <- merge(RECOV,DATA.4.one.state,by="date",all=TRUE)
}
q = which(RECOV[,"date"] == "2020-09-30")
RECOV = RECOV[1:q,]

# ===========================================
# Population:
#============================================
population2019 <- read.csv("Data_files/Population/population2019.csv", header=FALSE, sep=";")
# Reorder

population <- data.frame(state=vector.of.states,population2019=NaN)
vector.of.states
for(i in 1:nb.states){# state of residence
  State <- vector.of.states[i]
  indic.state <- which(State==population2019$V2)
  population$population2019[i] <- population2019$V3[indic.state]
}

# ===========================================
# Flow matrices:

# Commuting:
Matrix.commuting <- matrix(0,nb.states,nb.states)
Commuting <- read.csv("Data_files/TravelMatrices/Commuting.csv", sep=";")
for(i in 1:nb.states){# state of residence
  State <- vector.of.states[i]
  for(j in 1:nb.states){# state of work
    if(j!=i){
      State.work <- vector.of.states[j]
      DATA.aux <- subset(Commuting,(StateCode_residence==State)&(StateCode_work==State.work))
      if(dim(DATA.aux)[2]>0){
        Matrix.commuting[i,j] <- sum(DATA.aux$Workers.in.Commuting.Flow,na.rm = TRUE)
      }
    }
  }
}

Matrix.commuting <- diag(1/population$population2019) %*% Matrix.commuting

# Migration

State_to_State_Migrations_Table_2018 <- read.csv("Data_files/TravelMatrices/State_to_State_Migrations_Table_2018.csv", sep=";")
list.of.states <- names(State_to_State_Migrations_Table_2018)

Matrix.migration <- matrix(0,nb.states,nb.states)
for(i in 1:nb.states){# state of residence
  State <- vector.of.states[i]
  indic.state.origin <- which(State == list.of.states)
  for(j in 1:nb.states){# state of work
    State.destination <- vector.of.states[j]
    indic.state.destination <- which(State.destination == list.of.states)
    if(j!=i){
      if(dim(DATA.aux)[2]>0){
        Matrix.migration[i,j] <- State_to_State_Migrations_Table_2018[indic.state.destination,indic.state.origin]
      }
    }
  }
}



# ===========================================
# Policy dummies:

# Travel bans:

travelban <- t(read.csv("Data_files/Policies/travelban.csv", sep=","))
travelban <- data.frame(travelban[2:nrow(travelban),])
dates.travelban <- rownames(travelban)[2:nrow(travelban)]
names.travelban <- travelban[1,]
names.travelban.corrected <- rep(NA, length(names.travelban))
for(i in 1:length(names.travelban)){
  names.travelban.corrected[i] <- as.character(names.travelban[1,i])
}
for(i in 1:length(dates.travelban)){
  dates.travelban[i] <- substr(dates.travelban[i],2, 11)
}
dates.travelban     <- as.Date(dates.travelban, format = "%Y.%m.%d")
travelban           <- travelban[2:nrow(travelban),]
rownames(travelban) <- dates.travelban
names(travelban)    <- names.travelban.corrected

DATA.travelban <- data.frame(dates = dates.travelban)
for(i in 1:nb.states){
  State                     <- vector.of.states[i]
  indic.state.in.travelban  <- which(names(travelban)==State)
  data.to.fill              <- as.matrix(travelban[,indic.state.in.travelban])
  data.to.fill[data.to.fill!=0] <- 1
  eval(parse(text = gsub(" ","",paste("DATA.travelban$indic.travelban.",
                                      State," <- as.numeric(as.matrix(data.to.fill))", sep=""))))
}



# Stay at home:

stayathome <- t(read.csv("Data_files/Policies/stayathome.csv", sep=","))
stayathome <- data.frame(stayathome[2:nrow(stayathome),])
dates.stayathome <- rownames(stayathome)[2:nrow(stayathome)]
names.stayathome <- stayathome[1,]
names.stayathome.corrected <- rep(NA, length(names.stayathome))
for(i in 1:length(names.stayathome)){
  names.stayathome.corrected[i] <- as.character(names.stayathome[1,i])
}
for(i in 1:length(dates.stayathome)){
  dates.stayathome[i] <- substr(dates.stayathome[i],2, 11)
}
dates.stayathome     <- as.Date(dates.stayathome, format = "%Y.%m.%d")
stayathome           <- stayathome[2:nrow(stayathome),]
rownames(stayathome) <- dates.stayathome
names(stayathome)    <- names.stayathome.corrected

DATA.stayathome <- data.frame(dates = dates.stayathome)
for(i in 1:nb.states){
  State <- vector.of.states[i]
  indic.state.in.stayathome <- which(names(stayathome)==State)
  eval(parse(text = gsub(" ","",paste("DATA.stayathome$indic.stayathome.",
                                      State," <- as.numeric(as.matrix(stayathome[,indic.state.in.stayathome]))", sep=""))))
}

# Loading the travel data from Gustavo 
#=================================================
interstate_travel_taf         <- read.csv("Data_files/TravelMatrices/interstate_travel_taf.csv")
percentage_outofstate_travel  <- read.csv("Data_files/TravelMatrices/percentage_outofstate_travel.csv", sep=";")



DATA.travel <- matrix(0, nb.states, nb.states)
where.in.data.is.state <- NULL
vector.state.prop.travelout <- NULL
for(i in 1:nb.states){
  State                   <- vector.of.states[i]
  
  state.prop.travelout    <- percentage_outofstate_travel[which(percentage_outofstate_travel[,1]==State), 2]
  where.in.data.is.state  <- c(where.in.data.is.state, which(substr(as.character(interstate_travel_taf[,1]), 6, 7) == State))
  state.direction.to      <- interstate_travel_taf[which(substr(as.character(interstate_travel_taf[,1]), 6, 7) == State),]
  
  DATA.travel[i,] <- state.prop.travelout * as.numeric(state.direction.to[2:(length(vector.of.states)+1)])
  
  vector.state.prop.travelout <- c(vector.state.prop.travelout,state.prop.travelout)
}

DATA.travel <- DATA.travel[,where.in.data.is.state]


# Stay at home
stayathome_reduction_mobility <- read.csv("Data_files/TravelMatrices/stayathome_reduction_mobility.csv")
DATA.reduction_mobility <- rep(NA, nb.states)
for(i in 1:nb.states){
  State                       <- vector.of.states[i]
  DATA.reduction_mobility[i]  <- stayathome_reduction_mobility[which(stayathome_reduction_mobility[,1]== State),2]
}



# Reduction mobility travelban
travelban_reduction_travel <- read.csv("Data_files/TravelMatrices/travelban_reduction_travel.csv")
DATA.reduction_travel <- rep(NA, nb.states)
for(i in 1:nb.states){
  State                       <- vector.of.states[i]
  DATA.reduction_travel[i]  <- travelban_reduction_travel[which(travelban_reduction_travel[,1]== State),2]
}


# Masks mandate
masks <- t(read.csv("Data_files/Policies/mask.csv"))
masks <- data.frame(masks[2:nrow(masks),])
dates.masks <- rownames(masks)[2:nrow(masks)]
names.masks <- masks[1,]
names.masks.corrected <- rep(NA, length(names.masks))
for(i in 1:length(names.masks)){
  names.masks.corrected[i] <- as.character(names.masks[1,i])
}
for(i in 1:length(dates.masks)){
  dates.masks[i] <- substr(dates.masks[i],2, 11)
}
dates.masks     <- as.Date(dates.masks, format = "%Y.%m.%d")
masks           <- masks[2:nrow(masks),]
rownames(masks) <- dates.masks
names(masks)    <- names.masks.corrected

DATA.masks <- data.frame(dates = dates.masks)
for(i in 1:nb.states){
  State <- vector.of.states[i]
  indic.state.in.masks <- which(names(masks)==State)
  eval(parse(text = gsub(" ","",paste("DATA.masks$indic.masks.",
                                      State," <- as.numeric(as.matrix(masks[,indic.state.in.masks]))", sep=""))))
}





vector.of.states.full <- matrix(c(
  'DISTRICT OF COLUMBIA', "DC",
  "ALABAMA",	"AL",
  'ALASKA',	"AK",
  'ARIZONA',	'AZ',
  'ARKANSAS',	'AR',
  'CALIFORNIA',	'CA',
  'COLORADO',	'CO',
  'CONNECTICUT',	'CT',
  'DELAWARE',	'DE',
  'FLORIDA',	'FL',
  'GEORGIA',	'GA',
  'HAWAII',	'HI',
  'IDAHO',	'ID',
  'ILLINOIS',	'IL',
  'INDIANA',	'IN',
  'IOWA',	'IA',
  'KANSAS',	'KS',
  'KENTUCKY',	'KY',
  'LOUISIANA',	'LA',
  'MAINE',	'ME',
  'MARYLAND',	'MD',
  'MASSACHUSETTS',	'MA',
  'MICHIGAN',	'MI',
  'MINNESOTA',	'MN',
  'MISSISSIPPI',	'MS',
  'MISSOURI',	'MO',
  'MONTANA',	'MT',
  'NEBRASKA',	'NE',
  'NEVADA',	'NV',
  'NEW HAMPSHIRE',	'NH',
  'NEW JERSEY',	'NJ',
  'NEW MEXICO',	'NM',
  'NEW YORK',	'NY',
  'NORTH CAROLINA',	'NC',
  'NORTH DAKOTA',	'ND',
  'OHIO',	'OH',
  'OKLAHOMA',	'OK',
  'OREGON',	'OR',
  'PENNSYLVANIA',	'PA',
  'RHODE ISLAND',	'RI',
  'SOUTH CAROLINA',	'SC',
  'SOUTH DAKOTA',	'SD',
  'TENNESSEE',	'TN',
  'TEXAS',	'TX',
  'UTAH',	'UT',
  'VERMONT',	'VT',
  'VIRGINIA',	'VA',
  'WASHINGTON',	'WA',
  'WEST VIRGINIA',	'WV',
  'WISCONSIN',	'WI',
  'WYOMING',	'WY'), ncol = 2, byrow = T)

ordered.state.names <- vector.of.states.full[match(vector.of.states, vector.of.states.full[,2]),1]
