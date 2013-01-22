function base_path, type

  
  if n_elements(type) ne 0 then begin
     case type of
        'data': path = '/Users/bryna/Documents/Physics/idl_data_files/'
        'plots': path = '/Users/bryna/Documents/Physics/idl_plots/'
        else: message, 'type modifier not recognized.'
     endcase
  endif else path = '/Users/bryna/Documents/Physics/hazelton_git/idl_stuff/'

  return, path
end
