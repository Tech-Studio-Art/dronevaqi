const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch();
  const p = await b.newPage();
  await p.goto('file://' + process.argv[2], { waitUntil: 'load' });
  await p.pdf({ path: process.argv[3], format: 'A4', printBackground: true,
    displayHeaderFooter: true, headerTemplate: '<span></span>',
    footerTemplate: '<div style="font-size:7pt;width:100%;text-align:center;color:#666">Drone V AQI pilot proposal · page <span class="pageNumber"></span> of <span class="totalPages"></span></div>',
    margin: { top: '15mm', bottom: '16mm', left: '15mm', right: '15mm' } });
  await b.close();
})();
