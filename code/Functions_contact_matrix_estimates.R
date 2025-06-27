#' 
#' Author: Imelda Trejo Lorenzo.  
#' Date: June 27, 2025
#' Project grant name: CCM-UNAM & SECIHTI, PEE-2025-G-293
#' 
#'     
#' Aggregate Contact Matrix by Age Subgroups FUNCTION is base on the formulas 
#' given in the appendix of the paper: 
#' Promoting healthy populations as a pandemic preparedness strategy: a simulation study from Mexico
#' https://doi.org/10.1016/j.lana.2024.100682
#'
#' Aggregates a fine-scale contact matrix (e.g., 16x16) into a reduced matrix (e.g., 4x4),
#' based on custom-defined age subgroups and weighted by population distribution.
#' 
#' The 16x16 contact matrix row/column names were adapted from the paper: Projecting contact matrices 
#' in 177 geographical regions by Prem K et. al. PLoS Comput Biol 20(9): e1012454. 
#' https://doi.org/10.1371/journal.pcbi.1012454  
#'
#' @param contact_matrix A data frame with contact data. Must include columns:
#'   - `location_contact`: name of the setting (e.g., "home", "work", "school", etc.)
#'   - `setting`: the overall context (this function uses rows where setting == "overall")
#'   - `age_contactor`: fine-scale age group of the contactor
#'   - `age_contactee`: fine-scale age group of the contactee
#'   - `mean_number_of_contacts`: mean number of contacts per person
#'
#' @param population A data frame with population sizes by fine-scale age group. 
#'   Must include:
#'   - `AGE_GROUP`: fine-scale age group
#'   - `POPULATION`: number of individuals in that age group
#'
#' @param subgroups A named list where each name is a subgroup label (e.g., "0-19") 
#'   and each element is a character vector of fine-scale age groups included in that subgroup.
#'
#' @return A long-format data frame with columns:
#'   - `location_contact`: setting name
#'   - `age_contactee`: aggregated contactee age group
#'   - `age_contactor`: aggregated contactor age group
#'   - `mean_number_of_contacts`: weighted mean number of contacts from contactor to contactee
#'
#' @examples
#' # Example usage:
#' aggregate_contact_matrix(contact_matrix, population, subgroups)
#' 
aggregate_contact_matrix <- function(contact_matrix, population, subgroups) {
  
  new_contact_matrix <- list()
  row_index <- 1
  
  for (location_name in unique(contact_matrix$location_contact)) {
    
    Matriz_location <- contact_matrix %>%
      filter(setting == "overall", location_contact == location_name) %>%
      select(age_contactee, age_contactor, mean_number_of_contacts, location_contact)
    
    for (group_name_contactor in names(subgroups)) {
      
      subgroup_contactor <- subgroups[[group_name_contactor]]
      
      pop_contactor <- population %>%
        filter(AGE_GROUP %in% subgroup_contactor)
      
      total_pop_contactor <- sum(pop_contactor$POPULATION)
      
      for (group_name_contactee in names(subgroups)) {
        
        subgroup_contactee <- subgroups[[group_name_contactee]]
        contact_sum <- 0
        
        for (age_group in subgroup_contactor) {
          
          age_pop <- pop_contactor %>%
            filter(AGE_GROUP == age_group) %>%
            pull(POPULATION)
          
          Matriz_sub <- Matriz_location %>%
            filter(age_contactor == age_group,
                   age_contactee %in% subgroup_contactee)
          
          contact_sum <- contact_sum + sum(Matriz_sub$mean_number_of_contacts) * age_pop
        }
        
        contact_mean <- contact_sum / total_pop_contactor
        
        new_contact_matrix[[row_index]] <- tibble(
          location_contact = location_name,
          age_contactee = group_name_contactee,
          age_contactor = group_name_contactor,
          mean_number_of_contacts = contact_mean
        )
        
        row_index <- row_index + 1
      }
    }
  }
  return(bind_rows(new_contact_matrix))
}


grouping_age<- function(age_breaks,age_labels,population){

#stratify the population by groups and return their population size
  
population_grouped <-population %>%
  group_by(EDAD)%>%
  mutate(
    AGE_GROUP = cut(
      EDAD,
      breaks = age_breaks,
      right = FALSE,
      labels = age_labels,
      include.lowest = TRUE
    )
  )%>%
  group_by(AGE_GROUP) %>%
  summarise(POPULATION = sum(POBLACION, na.rm = TRUE)) 
}


contact_matrix_location <- function(df_data,location_name){
 #return a matrix four by four in the location_name place 
  
  names_subgroups <- unique(df_data$age_contactee)
  
  matrix <- matrix(0, nrow = length(names_subgroups), ncol = length(names_subgroups))
  
  # Índices inician en 1
  index_row <- 1
  
  for (group_name_contactor in names_subgroups) {
    
    index_col <- 1  # Reiniciar por cada fila
    
    for (group_name_contactee in names_subgroups) {
      
      Matrix_df <- df_data %>%
        filter(location_contact == location_name,
               age_contactee == group_name_contactee,
               age_contactor == group_name_contactor)
      
      matrix[index_row, index_col] <- Matrix_df$mean_number_of_contacts
      index_col <- index_col + 1
    }
    
    index_row <- index_row + 1
  }
  return(matrix)
}



#' Compute Consistent Symmetric Contact Matrix
#'
#'Algorithm base on Imelda Trejo's Paper
#'
#' Given an asymmetric contact matrix stratified by age and a corresponding population vector,
#' this function symmetrizes the matrix by averaging contacts from i→j and j→i, weighted by 
#' the size of the contactor population group.
#'
#' @param contact_matrix A data frame with contact information. Must include:
#'   - `location_contact`: name of the setting (e.g., "home", "work", "all", etc.)
#'   - `age_contactor`: age group initiating contact
#'   - `age_contactee`: age group receiving contact
#'   - `mean_number_of_contacts`: average number of contacts per person
#'
#' @param population A data frame with population sizes by age group.
#'   Must include:
#'   - `AGE_GROUP`: age group (matching `age_contactor` and `age_contactee`)
#'   - `POPULATION`: total population in that group
#'
#' @return A long-format data frame with columns:
#'   - `location_contact`, `age_contactee`, `age_contactor`, `mean_number_of_contacts` (symmetric)
#'
#' @examples
#' compute_consistent_matrix(contact_matrix, population)
#' 
#' 
compute_consistent_matrix <- function(contact_matrix, population) {
  
  consistent_contact_matrix <- list()
  row_index <- 1
  names_subgroups <- unique(contact_matrix$age_contactee)
  
  for (location_name in unique(contact_matrix$location_contact)) {
    
    for (group_name_contactor in names_subgroups) {
      
      pop_contactor <- population %>%
        filter(AGE_GROUP == group_name_contactor) %>%
        pull(POPULATION)
      
      for (group_name_contactee in names_subgroups) {
        
        pop_contactee <- population %>%
          filter(AGE_GROUP == group_name_contactee) %>%
          pull(POPULATION)
        
        contactor_to_contactee <- contact_matrix %>%
          filter(location_contact == location_name,
                 age_contactee == group_name_contactee,
                 age_contactor == group_name_contactor) %>%
          pull(mean_number_of_contacts)
        
        contactee_to_contactor <- contact_matrix %>%
          filter(location_contact == location_name,
                 age_contactee == group_name_contactor,
                 age_contactor == group_name_contactee) %>%
          pull(mean_number_of_contacts)
        
        # Compute consistent average (symmetric, population-weighted)
        consistent_contact_mean <- (
          contactor_to_contactee * pop_contactor +
            contactee_to_contactor * pop_contactee
        ) / (2 * pop_contactor)
        
        consistent_contact_matrix[[row_index]] <- tibble(
          location_contact = location_name,
          age_contactee = group_name_contactee,
          age_contactor = group_name_contactor,
          mean_number_of_contacts = consistent_contact_mean
        )
        
        row_index <- row_index + 1
      }
    }
  }
  
  bind_rows(consistent_contact_matrix)
}


