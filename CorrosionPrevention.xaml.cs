using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices.WindowsRuntime;
using System.Threading.Tasks;
using System.Xml.Linq;
using Windows.Foundation;
using Windows.Foundation.Collections;
using Windows.Security.Cryptography;
using Windows.Storage;
using Windows.Storage.Streams;
using Windows.UI.Popups;
using Windows.UI.Xaml;
using Windows.UI.Xaml.Controls;
using Windows.UI.Xaml.Controls.Primitives;
using Windows.UI.Xaml.Data;
using Windows.UI.Xaml.Documents;
using Windows.UI.Xaml.Input;
using Windows.UI.Xaml.Media;
using Windows.UI.Xaml.Navigation;

// The Blank Page item template is documented at https://go.microsoft.com/fwlink/?LinkId=234238

namespace CorrosionModel
{

    /// <summary>
    /// An empty page that can be used on its own or navigated to within a Frame.
    /// </summary>
    public sealed partial class CorrosionPrevention : Page
    {
        public StorageFile treatmentFile;
        public string AnodeCategory = "Z", CathodeCategory = "Z";
        public string AnodeOptions = null, CathodeOptions = null, CorrosionPreventionOptions = null;

        public CorrosionPrevention()
        {
            this.InitializeComponent();
            
        }

        private void ReturnButton_Click(object sender, RoutedEventArgs e)
        {
            this.Frame.Navigate(typeof(MainPage));
        }

        protected override void OnNavigatedTo(NavigationEventArgs e)
        {
            CorrosionRecommendationPayload passedPayload = e.Parameter as CorrosionRecommendationPayload;
            treatmentFile = passedPayload.aFile;
            AnodeCategory = passedPayload.AnCat;
            CathodeCategory = passedPayload.CatCat;

            DisplayRecommendations();

            base.OnNavigatedTo(e);
        }

        private async void DisplayRecommendations()
        {
            if (treatmentFile != null && AnodeCategory != "Z" && CathodeCategory != "Z")
            {
                try
                {
                    using (IRandomAccessStream readstream = await treatmentFile.OpenAsync(FileAccessMode.Read))
                    {
                        ulong size64 = readstream.Size;
                        if (size64 <= uint.MaxValue)
                        {

                            uint size32 = (uint)size64;
                            IBuffer buffer = new Windows.Storage.Streams.Buffer(size32);
                            buffer = await readstream.ReadAsync(buffer, size32, InputStreamOptions.None);
                            string filecontent = GetStringFromBuffer(buffer);

                            XElement xmlroot = XElement.Parse(filecontent);
                            IEnumerable<XElement> TreatmentNames_Anode = from el in xmlroot.Elements("Data")
                                                                         where (string)el.Element("ActivityCategory").Value == AnodeCategory
                                                                         select el;

                            AnodeRecommendationTextBlock.TextIndent = 0;
                            AnodeRecommendationTextBlock.TextAlignment = TextAlignment.Justify;
                            AnodeRecommendationTextBlock.FontSize = 14;

                            foreach (XElement el in TreatmentNames_Anode)
                            {
                                IEnumerable<XElement> AvailOptions = from an_el in el.Elements("Item") select an_el;
                                foreach (XElement opt1 in AvailOptions)
                                {                                    
                                    Run run1 = new Run
                                    {
                                        Text = opt1.Value + Environment.NewLine// "/n/r"; // "&#x0d;&#x0a";
                                    };
                                    Paragraph paragraph1 = new Paragraph();
                                    paragraph1.Inlines.Add(run1);
                                    AnodeRecommendationTextBlock.Blocks.Add(paragraph1);
                                }
                            }
                                                     

                            IEnumerable<XElement> TreatmentNames_Cathode = from el in xmlroot.Elements("Data")
                                                                           where (string)el.Element("ActivityCategory").Value == CathodeCategory
                                                                           select el;
                            CathodeRecommendationTextBlock.TextIndent = 0;
                            CathodeRecommendationTextBlock.TextAlignment = TextAlignment.Justify;
                            CathodeRecommendationTextBlock.FontSize = 14;

                            foreach (XElement el in TreatmentNames_Cathode)
                            {
                                IEnumerable<XElement> AvailOptions = from an_el in el.Elements("Item") select an_el;
                                foreach (XElement opt1 in AvailOptions)
                                {
                                    //run3.Text += opt1.Value + Environment.NewLine;//"/n/r"; // "&#x0d;&#x0a";
                                    Run run3 = new Run
                                    {
                                        Text = opt1.Value + Environment.NewLine// + Environment.NewLine"/n/r"; // "&#x0d;&#x0a";
                                    };
                                    Paragraph paragraph3 = new Paragraph();
                                    paragraph3.Inlines.Add(run3);
                                    CathodeRecommendationTextBlock.Blocks.Add(paragraph3);
                                }
                            }
                            //paragraph3.Inlines.Add(run3);
                            //CathodeRecommendatioTextBlock.Blocks.Add(paragraph3);
                        }
                        else
                        {
                            await DisplayMessageAsync("Alert message", "The treatment data file is too large.  It needs to be < 4 GB.", "notification");
                        }
                    }
                }
                catch
                {
                    await DisplayMessageAsync("Fail message", "The treatment data file did not load correctly.", "notification");
                }
            }
        }

        static internal string GetStringFromBuffer(IBuffer buffer)
        {
            return CryptographicBuffer.ConvertBinaryToString(BinaryStringEncoding.Utf8, buffer);
        }

        public async Task DisplayMessageAsync(string title, string content, string dialogType)
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
