/* Exercise MoSt's real input readers without creating an X display. */
#define main most_program_main
#include "most.c"
#undef main

int main(int argc, char **argv)
{
   double *grid;
   int i;

   if (argc != 3) return 2;
   if (!strcmp(argv[1], "timing"))
   {
      struct SelStruct timing[2] = {{0}};
      int found = 0;
      InitNamelist();
      for (i = 0; i < NL_items; ++i)
      {
         if (NL_list[i].model == PLASIM &&
             (!strcmp(NL_list[i].name, "MPSTEP") ||
              !strcmp(NL_list[i].name, "NAOMOD")))
         {
            if (found >= 2) return 6;
            timing[found].Item = &NL_list[i];
            timing[found].type = SEL_INT;
            timing[found].iv = NL_list[i].ival;
            strcpy(timing[found].text, NL_list[i].name);
            found++;
         }
      }
      if (found != 2) return 6;
      ComEnd = &SelStart;
      SelStart.Next = &timing[0];
      timing[0].Next = &timing[1];
      Model = PLASIM;
      Planet = EARTH;
      Lsg = Oce = Cores = Multirun = 1;
      Resolution = atoi(argv[2]);
      if (Resolution < 1 || Resolution >= RESOLUTIONS) return 7;
      CheckPlasimNamelist();
      for (i = 0; i < 2; ++i) sprintf(timing[i].teva, "%d", timing[i].iv);
      WriteNamelistFile("plasim", 0);
      WriteNamelistFile("oceanmod", 0);
   }
   else if (!strcmp(argv[1], "settings"))
   {
      strcpy(SelStart.text, "TEST");
      SelStart.type = SEL_INT;
      SelStart.iv = -1;
      printf("READ=%d\n", ReadSettings(argv[2]));
      printf("VALUE=%d\n", SelStart.iv);
   }
   else if (!strcmp(argv[1], "defaults"))
   {
      InitSelections();
      printf("NOUTPUT=%d\n", noutput);
      printf("NGUI=%d\n", ngui);
      printf("NNETCDF=%d\n", nnetcdf);
      printf("NOUTPUT_AVG=%d\n", noutput_avg);
      printf("SEL_NETCDF_NO=%d\n", SelNetCDF ? SelNetCDF->no : -1);
      printf("SEL_NETCDF_AVG_NO=%d\n", SelNetCDFAvg ? SelNetCDFAvg->no : -1);
      ChangeModel(PLASIM);
      printf("PLASIM_SEL_NETCDF_NO=%d\n", SelNetCDF ? SelNetCDF->no : -1);
      printf("PLASIM_SEL_NETCDF_AVG_NO=%d\n", SelNetCDFAvg ? SelNetCDFAvg->no : -1);
      ChangeModel(PUMA);
      printf("PUMA_SEL_NETCDF_NO=%d\n", SelNetCDF ? SelNetCDF->no : -1);
      printf("PUMA_SEL_NETCDF_AVG_NO=%d\n", SelNetCDFAvg ? SelNetCDFAvg->no : -1);
   }
   else if (!strcmp(argv[1], "resolutions"))
   {
      struct SelStruct *Sel;
      InitSelections();
      ChangeModel(PLASIM);
      for (i = RES_T21, Sel = SelRes; i <= RES_T63; ++i, Sel = Sel->Next)
      {
         printf("RES_%d=%s:%d\n", i, Sel->text, Sel->iv);
      }
      ChangeResolution(RES_T63);
      printf("AFTER_CHANGE_RES=%d\n", Resolution);
      for (i = RES_T21, Sel = SelRes; i <= RES_T63; ++i, Sel = Sel->Next)
      {
         printf("AFTER_RES_%d=%s:%d\n", i, Sel->text, Sel->iv);
      }
      CheckPlasimNamelist();
      printf("LATITUDES=%d\n", Latitudes);
      ChangeModel(CAT);
      for (i = RES_T21, Sel = SelRes; i <= RES_T63; ++i, Sel = Sel->Next)
      {
         printf("CAT_RES_%d_NO=%d\n", i, Sel->no);
      }
   }
   else if (!strcmp(argv[1], "config_t63"))
   {
      InitSelections();
      ChangeModel(PLASIM);
      ReadSettings(argv[2]);
      UpdateResolution();
      CheckPlasimNamelist();
      printf("CFG_RESOLUTION=%d\n", Resolution);
      printf("CFG_LATITUDES=%d\n", Latitudes);
   }
   else
   {
      int latitudes = !strcmp(argv[1], "reference32") ? 32 :
                      !strcmp(argv[1], "reference64") ? 64 : 4;
      grid = ReadGrid(argv[2], latitudes, 2.0);
      if (!grid) return 3;
      for (i = 0; i < latitudes; ++i)
      {
         int j;
         int stride = 2*latitudes + 1;
         for (j = 0; j < 2*latitudes; ++j)
         {
            if (!isfinite(grid[j + i*stride])) return 4;
            if (latitudes == 4 && grid[j + i*stride] != 2.0 * (j + i*8 + 1)) return 4;
         }
         if (grid[i*stride + 2*latitudes] != grid[i*stride]) return 5;
      }
      free(grid);
   }
   return 0;
}
