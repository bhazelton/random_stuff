pro read_metafits, filename, tiledata = tiledata, flagged_tilepol = flagged_tilepol


  data=readfits(filename, header2, exten=1)
  
  tbinfo, header, tb_str
  
  input = tbget(tb_str, data, 'input')
  antenna = tbget(tb_str, data, 'antenna')
  tile = tbget(tb_str, data, 'tile')
  pol = tbget(tb_str, data, 'pol')
  rx = tbget(tb_str, data, 'rx')
  slot = tbget(tb_str, data, 'slot')
  flag = tbget(tb_str, data, 'flag')
  length = tbget(tb_str, data, 'length')
  north = tbget(tb_str, data, 'north')
  east = tbget(tb_str, data, 'east')
  height = tbget(tb_str, data, 'height')
  
  tiledata = {input:input, antenna:antenna, tile:tile, pol:pol, rx:rx, slot:slot, flag:flag, $
    length:length, north:north, east:east, height:height}
    
  wh_flag = where(flag gt 0, count_flag)
  if count_flag gt 0 then begin
    for i=0, count_flag-1 do if i eq 0 then flagged_tilepol = number_formatter(tile[wh_flag[i]]) + pol[wh_flag[i]] $
    else flagged_tilepol = [flagged_tilepol, number_formatter(tile[wh_flag[i]]) + pol[wh_flag[i]]]
  endif
  
  
  
end