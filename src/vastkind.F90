      module vast_kind_param                                        
         integer, parameter :: byte_log = selected_int_kind(2)      
         integer, parameter :: short_log = selected_int_kind(4)     
         integer, parameter :: long_log = selected_int_kind(18)     
         integer, parameter :: byte = selected_int_kind(2)          
         integer, parameter :: short = selected_int_kind(4)         
         integer, parameter :: long = selected_int_kind(18)         
         integer, parameter :: double = selected_real_kind(8)      
         integer, parameter :: extended = selected_real_kind(30)    
         integer, parameter :: double_ext = selected_real_kind(50)  
         integer, parameter :: dble_complex = selected_real_kind(14)
         integer, parameter :: ext_complex = selected_real_kind(30) 
      end module vast_kind_param                                    
