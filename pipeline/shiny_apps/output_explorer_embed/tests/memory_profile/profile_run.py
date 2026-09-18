"""Drive one profile run: upload bundle -> walk every tab -> screenshot each -> quit.

Assumes the Shiny app is already running on the given port and a sampler is
running against that R process. Writes per-tab checkpoint timestamps to a
newline-delimited JSON file so the sampler CSV can be annotated afterwards.

Firefox + geckodriver locations are read from FIREFOX_BIN / GECKODRIVER_BIN
env vars if set, otherwise the ones on PATH are used.

Usage:
  python profile_run.py \
    --url http://127.0.0.1:4321 \
    --bundle path/to/bundle.tar.gz \
    --outdir path/to/results/N04
"""
import argparse, os, json, time, shutil
from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

ap = argparse.ArgumentParser()
ap.add_argument("--url", required=True)
ap.add_argument("--bundle", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--tab-wait", type=float, default=5.0, help="seconds to wait after clicking each tab")
args = ap.parse_args()

FIREFOX_BIN     = os.environ.get("FIREFOX_BIN")     or shutil.which("firefox")
GECKODRIVER_BIN = os.environ.get("GECKODRIVER_BIN") or shutil.which("geckodriver")
if not FIREFOX_BIN or not GECKODRIVER_BIN:
    raise SystemExit("firefox / geckodriver not found; set FIREFOX_BIN and GECKODRIVER_BIN env vars")

os.makedirs(args.outdir, exist_ok=True)
checkpoints_path = os.path.join(args.outdir, "checkpoints.jsonl")
cp_f = open(checkpoints_path, "w")

def cp(label, extra=None):
    rec = {"epoch_ms": int(time.time() * 1000), "label": label}
    if extra: rec.update(extra)
    cp_f.write(json.dumps(rec) + "\n"); cp_f.flush()
    print(f"[{time.strftime('%H:%M:%S')}] {label}" + (f" {extra}" if extra else ""))

opts = Options()
opts.add_argument("--headless")
opts.add_argument("--width=1440")
opts.add_argument("--height=1200")
opts.binary_location = FIREFOX_BIN
service = Service(executable_path=GECKODRIVER_BIN)
driver = webdriver.Firefox(options=opts, service=service)
driver.set_window_size(1440, 1200)

def shot(name):
    p = os.path.join(args.outdir, f"{name}.png")
    driver.save_screenshot(p)
    return p

def click_tab(title):
    # Use JS click so we bypass Selenium's "scroll into view" precheck, which
    # fails when the tab list overflows / when a modal alert is queued.
    driver.execute_script("""
      var t = arguments[0];
      var links = document.querySelectorAll('#main_tabs > li > a');
      for (var i=0; i<links.length; i++) {
        if (links[i].textContent.trim() === t) { links[i].click(); return true; }
      }
      throw new Error('tab not found: ' + t);
    """, title)

def dismiss_alerts():
    try:
        while True:
            a = driver.switch_to.alert
            a.accept()
    except Exception:
        pass

TABS = [
    "Overview",
    "GWAS QC",
    "SNP-h² & rG",
    "SNP Associations",
    "Molecular Associations",
    "Enrichment Analysis",
    "References",
    "Configuration",
]

try:
    cp("start")
    driver.get(args.url)
    WebDriverWait(driver, 30).until(
        EC.presence_of_element_located((By.CSS_SELECTOR, "input[type='file']"))
    )
    cp("page_loaded")
    shot("00_page_loaded")

    # Upload
    fi = driver.find_element(By.CSS_SELECTOR, "input[type='file']")
    fi.send_keys(args.bundle)
    cp("upload_sent", {"bundle": args.bundle, "size_bytes": os.path.getsize(args.bundle)})

    # Wait for app-ready flag set by config_flags observer (see app.R)
    WebDriverWait(driver, 300).until(
        lambda d: d.execute_script(
            "return document.documentElement.getAttribute('data-app-ready') === '1';"
        )
    )
    cp("app_ready")
    shot("01_app_ready")
    time.sleep(3)

    for i, tab in enumerate(TABS, start=1):
        try:
            dismiss_alerts()
            click_tab(tab)
            time.sleep(args.tab_wait)
            dismiss_alerts()
            cp(f"tab_{i:02d}_{tab.replace(' ', '_').replace('&','and')}")
            shot(f"tab_{i:02d}_{tab.replace(' ', '_').replace('&','and')}")
        except Exception as e:
            cp(f"tab_error_{tab}", {"err": str(e)[:200]})

    def walk_subtabs(section_name, css_selector):
        try:
            dismiss_alerts()
            click_tab(section_name)
            time.sleep(2)
            dismiss_alerts()
            titles = driver.execute_script("""
              var sel = arguments[0];
              return Array.from(document.querySelectorAll(sel)).map(function(a){return a.textContent.trim();}).filter(function(t){return t.length>0;});
            """, css_selector)
            seen = set()
            for t in titles:
                if t in seen: continue
                seen.add(t)
                try:
                    driver.execute_script("""
                      var sel = arguments[0]; var t = arguments[1];
                      var links = document.querySelectorAll(sel);
                      for (var i=0; i<links.length; i++) {
                        if (links[i].textContent.trim() === t) { links[i].click(); return; }
                      }
                    """, css_selector, t)
                    time.sleep(3)
                    dismiss_alerts()
                    cp(f"{section_name[:3].lower()}_sub_{t}")
                except Exception as e:
                    cp(f"{section_name[:3].lower()}_sub_error_{t}", {"err": str(e)[:150]})
            shot(f"99_{section_name[:3].lower()}_final")
        except Exception as e:
            cp(f"{section_name[:3].lower()}_walk_error", {"err": str(e)[:200]})

    walk_subtabs("Molecular Associations", "div[id^='mol_assoc'] ul.nav a")
    walk_subtabs("Enrichment Analysis",    "div[id^='enrichment'] ul.nav a")

    cp("done")

finally:
    cp_f.close()
    driver.quit()
