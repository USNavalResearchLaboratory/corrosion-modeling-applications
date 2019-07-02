using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices.WindowsRuntime;
using System.Threading;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Linq;
using Windows.Foundation;
using Windows.Foundation.Collections;
using Windows.UI.Xaml;
using Windows.UI.Xaml.Controls;
using Windows.UI.Xaml.Controls.Primitives;
using Windows.UI.Xaml.Data;
using Windows.UI.Xaml.Input;
using Windows.UI.Xaml.Media;
using Windows.UI.Xaml.Navigation;
using Windows.Storage;
using Windows.Storage.Streams;
using Windows.Security.Cryptography;
using Windows.UI.Popups;
using System.Diagnostics;
using Windows.UI;
using Windows.UI.Xaml.Documents;

// The Blank Page item template is documented at https://go.microsoft.com/fwlink/?LinkId=402352&clcid=0x409

namespace CorrosionModel
{
    public class CorrosionRecommendationPayload
    {
        public string AnCat { get; set; }
        public string CatCat { get; set; }
        public StorageFile aFile { get; set; }
    }

    /// <summary>
    /// An empty page that can be used on its own or navigated to within a Frame.
    /// </summary> 
    public sealed partial class MainPage : Page
    {
        public double AnodePotentialValue, CathodePotentialValue;
        public StorageFile sampleFile, treatmentFile;  
        public string[,] CompatibilityTable = new string[20, 61];
        public string AnodeCategory = "Z", CathodeCategory = "Z";

        //public string AnodeOptions = null, CathodeOptions = null, CorrosionPreventionOptions = null;

        public class Material
        {
            public string MaterialName { get; set; }
        }

        public MainPage()
        {
            this.InitializeComponent();
            LoadCompatibilityFile();
        }
               
        private async void CalculateButton_Click(object sender, RoutedEventArgs e)
        {
            string AnMat, CatMat;
            
            bool IndustrialAtmosphereCompatibility = false, MarineAtmosphereCompatibility = false, SeawaterCompatibility = false;

            string anode_val, cathode_val; 

            AnMat = Convert.ToString(ListOfAnodes.SelectionBoxItem);
            CatMat = Convert.ToString(ListOfCathodes.SelectionBoxItem);

            double an_val = 0.0, cat_val = 0.0, diff_val = 0.0;

            if (sampleFile != null)
            {
                try
                {
                    using (IRandomAccessStream readstream = await sampleFile.OpenAsync(FileAccessMode.Read))
                    {
                        ulong size64 = readstream.Size;
                        if (size64 <= uint.MaxValue)
                        {

                            uint size32 = (uint)size64;
                            IBuffer buffer = new Windows.Storage.Streams.Buffer(size32);
                            buffer = await readstream.ReadAsync(buffer, size32, InputStreamOptions.None);
                            string filecontent = GetStringFromBuffer(buffer);

                            XElement xmlroot = XElement.Parse(filecontent);
                            IEnumerable<XElement> check_for_anode_name = from el in xmlroot.Elements("Data")
                                                          where (string)el.Element("Name") == AnMat
                                                          select el;
                            foreach (XElement el in check_for_anode_name)
                            {
                                anode_val = el.Element("PotentialValue").Value;
                                AnodeCategory = el.Element("ActivityCategory").Value;
                                an_val = Convert.ToDouble(anode_val);
                                AnodePotential.Text = anode_val;
                            }

                            IEnumerable<XElement> check_for_cathode_name = from el in xmlroot.Elements("Data")
                                                                         where (string)el.Element("Name") == CatMat
                                                                           select el;
                            foreach (XElement el in check_for_cathode_name)
                            {
                                cathode_val = el.Element("PotentialValue").Value;
                                CathodeCategory = el.Element("ActivityCategory").Value;
                                cat_val = Convert.ToDouble(cathode_val);
                                CathodePotential.Text = cathode_val;
                            }

                            diff_val = cat_val - an_val;
                            PotentialDifference.Text = diff_val.ToString();

                            
                            int cat_search_val = 0;
                            if (cat_val > an_val)
                            {
                                cat_search_val = ReturnIntFromLetter(CathodeCategory);

                                for (int i = 0; i <= CompatibilityTable.GetLength(0)-1;i++)
                                {
                                    if (CompatibilityTable[i,0] == AnodeCategory)
                                    {
                                        int start_val = (cat_search_val * 3) + 1;
                                        if (CompatibilityTable[i, start_val]  == "C")
                                        {
                                            IndustrialAtmosphereCompatibility = true;
                                        }
                                        if (CompatibilityTable[i, start_val+1] == "C")
                                        {
                                            SeawaterCompatibility = true;
                                        }
                                        if (CompatibilityTable[i, start_val+2] == "C")
                                        {
                                            MarineAtmosphereCompatibility = true;
                                        }
                                    }
                                }
                            }
                            else
                            {
                                cat_search_val = ReturnIntFromLetter(AnodeCategory);
                                for (int i = 0; i <= CompatibilityTable.GetLength(0) - 1; i++)
                                {
                                    if (CompatibilityTable[i, 0] == CathodeCategory)
                                    {
                                        int start_val = (cat_search_val * 3) + 1;
                                        if (CompatibilityTable[i, start_val] == "C")
                                        {
                                            IndustrialAtmosphereCompatibility = true;
                                        }
                                        if (CompatibilityTable[i, start_val + 1] == "C")
                                        {
                                            SeawaterCompatibility = true;
                                        }
                                        if (CompatibilityTable[i, start_val + 2] == "C")
                                        {
                                            MarineAtmosphereCompatibility = true;
                                        }
                                    }
                                }
                            }

                            SolidColorBrush IASolidColorBrush = new SolidColorBrush();
                            SolidColorBrush SWSolidColorBrush = new SolidColorBrush();
                            SolidColorBrush MASolidColorBrush = new SolidColorBrush();


                            if (AnMat == CatMat)
                            {
                                IASolidColorBrush.Color = Color.FromArgb(255, 0, 255, 0);
                                IA_Rectangle.Fill = IASolidColorBrush;
                                SW_Rectangle.Fill = IASolidColorBrush;
                                MA_Rectangle.Fill = IASolidColorBrush;
                            }
                            else
                            {
                                if (IndustrialAtmosphereCompatibility)
                                {
                                    IASolidColorBrush.Color = Color.FromArgb(255, 0, 255, 0);
                                    IA_Rectangle.Fill = IASolidColorBrush;
                                }
                                else
                                {
                                    IASolidColorBrush.Color = Color.FromArgb(255, 255, 0, 0);
                                    IA_Rectangle.Fill = IASolidColorBrush;
                                }
                                if (SeawaterCompatibility)
                                {
                                    SWSolidColorBrush.Color = Color.FromArgb(255, 0, 255, 0);
                                    SW_Rectangle.Fill = SWSolidColorBrush;
                                }
                                else
                                {
                                    SWSolidColorBrush.Color = Color.FromArgb(255, 255, 0, 0);
                                    SW_Rectangle.Fill = SWSolidColorBrush;
                                }
                                if (MarineAtmosphereCompatibility)
                                {
                                    MASolidColorBrush.Color = Color.FromArgb(255, 0, 255, 0);
                                    MA_Rectangle.Fill = MASolidColorBrush;
                                }
                                else
                                {
                                    MASolidColorBrush.Color = Color.FromArgb(255, 255, 0, 0);
                                    MA_Rectangle.Fill = MASolidColorBrush;
                                }
                                
                            }
                                                       
                        }
                        else
                        {
                            await displayMessageAsync("Alert message", "Something went wrong with loading the potential data.", "notification");
                        }
                    }

                    CorrosionPreventionButton.IsEnabled = true;
                }
                catch
                {
                    await displayMessageAsync("Fail message", "The default data file did not load correctly.", "notification");
                }
            }
             
        }

        private void ListOfAnodes_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            ListOfCathodes.IsEnabled = true;
        }

        private void ListOfCathodes_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            CalculateButton.IsEnabled = true;           
        }

        private async void AppBar_LoadDefault_Button_Click(object sender, RoutedEventArgs e)
        {
            string currentDirectory = Directory.GetCurrentDirectory(); //$"C:\\Users\\steve\\Documents"; // 

            // Load the XML file from our project directory containing the purchase orders
            string potentialdata_filename = "SeawaterPotentialData.xml";          
            string galvPotentialFilepath = Path.Combine(currentDirectory, potentialdata_filename);
            StorageFile temp_file = await StorageFile.GetFileFromPathAsync(galvPotentialFilepath);

            string alloytreatment_filename = "AlloyTreatmentRecommendations.xml";
            string treatmentFilepath = Path.Combine(currentDirectory, alloytreatment_filename);
            StorageFile temp_file_2 = await StorageFile.GetFileFromPathAsync(treatmentFilepath);

            if (temp_file != null) { sampleFile = temp_file; }
            if (sampleFile != null)
            {
                try
                {
                    using (IRandomAccessStream readstream = await sampleFile.OpenAsync(FileAccessMode.Read))
                    {
                        ulong size64 = readstream.Size;
                        if (size64 <= uint.MaxValue)
                        {

                            uint size32 = (uint)size64;
                            IBuffer buffer = new Windows.Storage.Streams.Buffer(size32);
                            buffer = await readstream.ReadAsync(buffer, size32, InputStreamOptions.None);
                            string filecontent = GetStringFromBuffer(buffer);

                            XElement xmlroot = XElement.Parse(filecontent);
                            IEnumerable<XElement> MaterialNames = from el in xmlroot.Elements("Data") select el;

                            foreach (XElement el in MaterialNames)
                            {
                                ListOfAnodes.Items.Add(el.Element("Name").Value);
                                ListOfCathodes.Items.Add(el.Element("Name").Value);
                            }

                            ListOfAnodes.IsEnabled = true;
                        }
                        else
                        {
                            await displayMessageAsync("Alert message", "Your data file is too large.  It needs to be < 4 GB.", "notification");
                        }
                    }
                }
                catch
                {
                    await displayMessageAsync("Fail message", "The default data file did not load correctly.", "notification");
                }
            }

            if (temp_file_2 != null) { treatmentFile = temp_file_2; }
        }

        private async void LoadCompatibilityFile()
        {
            string filename = "TableI_Data.csv";
            string currentDirectory = Directory.GetCurrentDirectory(); //$"C:\\Users\\steve\\Documents"; // 
            string CompatibilityFilepath = Path.Combine(currentDirectory, filename);
            //StorageFolder storageFolder = await StorageFolder.GetFolderFromPathAsync(currentDirectory);
            StorageFile temp_file = await StorageFile.GetFileFromPathAsync(CompatibilityFilepath);

            if (temp_file != null)
            {
                try
                {
                    using (IRandomAccessStream readstream = await temp_file.OpenAsync(FileAccessMode.Read))
                    {
                        ulong size64 = readstream.Size;
                        if (size64 <= uint.MaxValue)
                        {

                            uint size32 = (uint)size64;
                            IBuffer buffer = new Windows.Storage.Streams.Buffer(size32);
                            buffer = await readstream.ReadAsync(buffer, size32, InputStreamOptions.None);
                            string filecontent = GetStringFromBuffer(buffer);

                            string[] data_from_comp_file = filecontent.Split(new char[] {',', '\r', '\n' }, StringSplitOptions.None);

                            int start_of_table = 63;
                            int num_rows_of_table = CompatibilityTable.GetLength(0);
                            int num_cols_of_table = CompatibilityTable.GetLength(1);

                            int array_index = 0;
                            int start_index = start_of_table, end_index = start_of_table + num_cols_of_table;
                            for (int i = 0; i <= num_rows_of_table - 1; i++)
                            {
                                array_index = start_index;
                                for (int j = 0; j <= num_cols_of_table - 1; j++)
                                {
                                    CompatibilityTable[i, j] = data_from_comp_file[array_index];
                                    array_index++;
                                }
                                start_index = end_index + 1;
                                end_index = start_index + num_cols_of_table;

                            }

                        }
                    }
                }
                catch
                {
                    await displayMessageAsync("Fail message", "The compatibility data file did not load correctly.", "notification");
                }
            }
        }

        private void CorrosionPreventionButton_Click(object sender, RoutedEventArgs e)
        {
            CorrosionRecommendationPayload payload = new CorrosionRecommendationPayload
            {
                AnCat = AnodeCategory,
                CatCat = CathodeCategory,
                aFile = treatmentFile
            };

            this.Frame.Navigate(typeof(CorrosionPrevention), payload);
        }

        private int ReturnIntFromLetter(string A)
        {
            int val = 0;

            switch (A)
            {
                case "A":
                    val = 0;
                    break;
                case "B":
                    val = 1;
                    break;
                case "C":
                    val = 2;
                    break;
                case "D":
                    val = 3;
                    break;
                case "E":
                    val = 4;
                    break;
                case "F":
                    val = 5;
                    break;
                case "G":
                    val = 6;
                    break;
                case "H":
                    val = 7;
                    break;
                case "I":
                    val = 8;
                    break;
                case "J":
                    val = 9;
                    break;
                case "K":
                    val = 10;
                    break;
                case "L":
                    val = 11;
                    break;
                case "M":
                    val = 12;
                    break;
                case "N":
                    val = 13;
                    break;
                case "O":
                    val = 14;
                    break;
                case "P":
                    val = 15;
                    break;
                case "Q":
                    val = 16;
                    break;
                case "R":
                    val = 17;
                    break;
                case "S":
                    val = 18;
                    break;
                case "T":
                    val = 19;
                    break;
            }

            return val;
        }

        private async void AppBar_LoadCustom_Button_Click(object sender, RoutedEventArgs e)
        {
            Windows.Storage.Pickers.FileOpenPicker picker = new Windows.Storage.Pickers.FileOpenPicker
            {
                ViewMode = Windows.Storage.Pickers.PickerViewMode.List
            };
            picker.FileTypeFilter.Add(".xml");
            //picker.SuggestedStartLocation = Windows.Storage.Pickers.PickerLocationId.Unspecified;
            StorageFile file = await picker.PickSingleFileAsync();
            if (file != null)
            {
                // Application now has read/write access to the picked file
                sampleFile = file;
                ListOfAnodes.IsEnabled = true;
            }
            else
            {
               await displayMessageAsync("Fail message","Your custom data file did not load correctly.", "notification");
            }
        }

        static internal string GetStringFromBuffer(IBuffer buffer)
        {
            return CryptographicBuffer.ConvertBinaryToString(BinaryStringEncoding.Utf8, buffer);
        }

        private void CorrosionPlot_Click(object sender, RoutedEventArgs e)
        {
            this.Frame.Navigate(typeof(CorrosionPrevention));
        }

        public async Task displayMessageAsync(string title, string content, string dialogType)
        {
            var messageDialog = new MessageDialog(content, title);
            if (dialogType == "notification")
            {
                //Do nothing here.Display normal notification MessageDialog
            }
            else
            {
                //Dipplay questions-Yes or No- MessageDialog
                messageDialog.Commands.Add(new UICommand("Yes", null));
                messageDialog.Commands.Add(new UICommand("No", null));
                messageDialog.DefaultCommandIndex = 0;
            }

            messageDialog.CancelCommandIndex = 1;
            var cmdResult = await messageDialog.ShowAsync();
            if (cmdResult.Label == "Yes")
            {
                Debug.WriteLine("My Dialog answer label is:: " + cmdResult.Label);
            }
            else
            {
                Debug.WriteLine("My Dialog answer label is:: " + cmdResult.Label);
            }
        }
    }

}
