function base_path, data=data, plots=plots

  if n_elements(data) ne 0 and n_elements(plots) ne 0 then message, 'data and plots options cannot be used simultaneously'

  if n_elements(data) ne 0 then path = '/Users/bryna/Documents/Physics/idl_data_files/'
  if n_elements(plots) ne 0 then path = '/Users/bryna/Documents/Physics/idl_plots/'
  if n_elements(path) ne 0 then path = '/Users/bryna/Documents/Physics/bryna_git/idl_stuff/'

  return, path
end
