function base_path, type


  if n_elements(type) ne 0 then begin
     case type of
        'data': path = '/Users/bryna/Projects/Physics/data_files/'
        'plots': path = '/Users/bryna/Projects/Physics/idl_plots/'
        else: message, 'type modifier not recognized.'
     endcase
  endif else path = '/Users/bryna/Projects/Physics/random_stuff/idl_stuff/'

  return, path
end
