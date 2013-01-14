
pro print_django_cmd, table = table

  tables = ['beamformer', 'cable', 'bfconn', 'rxbfconn']
  if n_elements(table) eq 0 then table = tables[0]
  wh = where(tables eq table, count)
  if count eq 0 then message, 'table not recognized'

  tiles = [11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 31, 32, 33, 34, 35, 36, 37, 38, 41, 42, 43, 44, 45, 46, $
           47, 48, 51, 52, 53, 54, 55, 56, 57, 58, 61, 62, 63, 64, 65, 66, 67, 68, 71, 72, 73, 74, 75, 76, 77, 78, 81, 82, 83, 84, $
           85, 86, 87, 88, 91, 92, 93, 94, 95, 96, 97, 98, 101, 102, 103, 104, 105, 106, 107, 108, 111, 112, 113, 114, 115, 116, 117, $
           118, 121, 122, 123, 124, 125, 126, 127, 128, 131, 132, 133, 134, 135, 136, 137, 138, 141, 142, 143, 144,  145, 146, 147, $
           148, 151, 152,  153, 154, 155, 156, 157, 158,  161, 162, 163, 164, 165, 166, 167, 168]  
 
  lengths = [90, 90, 150, 150, 150, 90, 150, 90, 90, 90, 90, 150, 150, 150, 150, 90, 90, 90, 90, 150, 150, 150, 150, 150, 150, $
             150, 150, 150, 150, 90, 150, 90, 90, 400, 90, 320, 524, 230, 230, 90, 150, 150, 150, 150, 230, 150, 150, 150, 230, $
             230, 400, 230, 230, 230, 90, 150, 230, 150, 230, 230, 230, 230, 230, 90, 320, 230, 150, 230, 230, 230, 230, 230, $
             400, 400, 150, 524, 230, 150, 524, 524, 524, 400, 524, 400, 524, 524, 524, 524, 524, 524, 400, 320, 524, 524, 400, $
             524, 524, 24, 524, 400, 320, 400, 320, 524, 524, 524, 320, 524, 524, 320, 400, 400, 524, 230, 524, 524, 400, 400, $
             400, 400, 90, 230, 524, 524, 524, 400, 400, 320]

  case table of
     'beamformer': begin
        for i=0, n_elements(tiles) -1 do begin
           print, 'b = Beamformer(beamformernumber = ' + number_formatter(tiles[i]) + ', sn = "BF' + number_formatter(tiles[i]) + $
                  '", gainx = 11, gainy = 11)'
           print, 'b.save()'
        endfor
     end
     'cable': begin

        for i=0, n_elements(tiles) -1 do begin
           if lengths[i] lt 320 then begin
              type = 'RG-6'
              prefix = 'RG'
           endif else begin
              type = 'LMR400-75'
              prefix = 'LMR'
           endelse

           print, 'c = Cable(name = "' + prefix + number_formatter(tiles[i]) + '", type = "' + type + '", eleclength = ' + $
                  number_formatter(lengths[i]) + ', physlength = ' + number_formatter(lengths[i]) + $
                  ', commDate = "2012-06-12 00:00:00")'
           print, 'c.save()'
        endfor
     end
     'bfconn':begin
        for i=0, n_elements(tiles) -1 do begin
           print, 'b = Beamformer.objects.get(beamformernumber = ' + number_formatter(tiles[i]) + ')'
           print, 't = Tile.objects.get(tilenumber = ' + number_formatter(tiles[i]) + ')'
           print, 'c = BeamformerTileConnection(beamformer = b, tile = t)'
           print, 'c.save()'
        endfor
     end
     'rxbfconn': begin
        for i=0, n_elements(tiles) -1 do begin
           if lengths[i] lt 320 then prefix = 'RG'else  prefix = 'LMR'
           rec_num = (i / 8) + 1
           slot_num = (i mod 8) + 1
           
           print, 'b = Beamformer.objects.get(beamformernumber = ' + number_formatter(tiles[i]) + ')'
           print, 'c = Cable.objects.get(name = "' + prefix + number_formatter(tiles[i]) + '")'
           print, 'r = Receiver.objects.get(receivernumber = ' + number_formatter(rec_num) + ')'
           print, 's = ' + number_formatter(slot_num)
           print, 'a = ReceiverBeamformerConnection(beamformer = b, cable = c, slot = s, receiver = r)'
           print, 'a.save()'
        endfor
     end
  endcase

end
