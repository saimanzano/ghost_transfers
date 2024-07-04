#!/usr/bin/env Rscript
#  Copyright Saioa Manzano-Morales, Moisès Bernabeu & Toni Gabaldón <saioa.manzano@bsc.es>
#  
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.  

# Functions definition
get_donors <- function(x) {
  table <- x$table
  res <- c()
  for (i in unique(table$iter)) {
    st <- table[which(table$iter == i), ]
    shift_prop <- table(st[, 'shift']) / 2
    if (is.na(shift_prop['TRUE'])) {
      shift_prop <- 0
    } else {
      shift_prop <- shift_prop['TRUE']
    }
    res <- rbind(res, data.frame('donor1' = st[1, 'donor'],
                                 'donor2' = st[2, 'donor'],
                                 'shift_prop' = shift_prop,
                                 'shift' = st[, 'shift'][1],
                                 'older_th_feca' = st[, 'older_th_feca'][1],
                                 'iter' = i))
  }
  row.names(res) <- NULL
  return(res)
}


## Data import
cat('Loading data\n')
load('../outputs/shift_list.RData')
cat('Data loaded\n')

cat('Initiating donors analysis\n')
donors_shift <- do.call(rbind, lapply(shiftl, get_donors))

cat('Writing donors table\n')
write.table(donors_shift, file = '../outputs/donors_shift.tsv', sep='\t')
cat('Table written\n')
