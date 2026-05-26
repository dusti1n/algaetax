import sys
from pathlib import Path

from PyQt6.QtCore import QProcess, QTimer
from PyQt6.QtGui import QColor, QIcon, QPainter, QPen
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QFileDialog,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProxyStyle,
    QPushButton,
    QStyle,
    QVBoxLayout,
    QWidget,
)


ICON_PATH = Path("documentation/images/algaetax_icon.png")
DEFAULT_PR2_PATH = "database/pr2_version_5.1.1_taxonomy.xlsx"
DEFAULT_NCBI_EMAIL = "anonymous@ncbi.com"
DEFAULT_ALGAEBASE_URL = "https://api.algaebase.org/v1.3/species"


class DarkCheckboxStyle(QProxyStyle):
    """Custom dark checkbox indicator with an X marker for checked states."""

    def drawPrimitive(self, element, option, painter, widget=None):
        if element == QStyle.PrimitiveElement.PE_IndicatorCheckBox:
            rect = option.rect.adjusted(1, 1, -1, -1)

            painter.save()
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
            painter.setPen(QPen(QColor("#59706C"), 2))
            painter.setBrush(QColor("#2a2f2e"))
            painter.drawRect(rect)

            if option.state & QStyle.StateFlag.State_On:
                painter.setPen(QPen(QColor("#f2f4f3"), 2))
                x_rect = rect.adjusted(5, 5, -5, -5)
                painter.drawLine(x_rect.topLeft(), x_rect.bottomRight())
                painter.drawLine(x_rect.topRight(), x_rect.bottomLeft())

            painter.restore()
            return

        super().drawPrimitive(element, option, painter, widget)


class AlgaeTaxGUI(QWidget):
    """Graphical user interface for creating configs and running algaetax."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("algaetax: taxonomic data query GUI")
        self.setWindowIcon(QIcon(str(ICON_PATH)))
        self.setFixedSize(1140, 820)

        self.process = None
        self.process_status_file = None
        self.process_timer = QTimer(self)
        self.process_timer.setInterval(2000)
        self.process_timer.timeout.connect(self.check_process_status)

        self.build_ui()
        self.apply_style()

    def build_ui(self):
        """Build the complete two-column GUI layout."""
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(18, 14, 18, 14)
        main_layout.setSpacing(20)

        left_column = QVBoxLayout()
        right_column = QVBoxLayout()
        left_column.setSpacing(18)
        right_column.setSpacing(18)

        main_layout.addLayout(left_column, 1)
        main_layout.addLayout(right_column, 1)

        self.build_left_column(left_column)
        self.build_right_column(right_column)

    def build_left_column(self, left_column: QVBoxLayout):
        """Build the title, general settings, options, and database selector."""
        title_row = QHBoxLayout()

        title = QLabel("algaetax: taxonomy data query GUI")
        title.setObjectName("title")
        title.setWordWrap(False)

        close_btn = QPushButton("Close")
        close_btn.setFixedWidth(90)
        close_btn.clicked.connect(self.close)

        title_row.addWidget(title)
        title_row.addStretch()
        title_row.addWidget(close_btn)
        left_column.addLayout(title_row)

        subtitle = QLabel("See the README or config file for parameter details.")
        subtitle.setWordWrap(True)
        subtitle.setObjectName("infoHint")
        left_column.addWidget(subtitle)

        left_column.addWidget(self.create_general_box())
        left_column.addWidget(self.create_taxa_box())
        left_column.addWidget(self.create_options_box())
        left_column.addWidget(self.create_database_box())
        left_column.addStretch()

    def build_right_column(self, right_column: QVBoxLayout):
        """Build database path, API, and control sections."""
        right_column.addWidget(self.create_pr2_box())
        right_column.addWidget(self.create_ncbi_box())
        right_column.addWidget(self.create_algaebase_box())
        right_column.addWidget(self.create_controls_box())
        right_column.addStretch()

    def create_general_box(self) -> QGroupBox:
        box = self.section_box("General")
        layout = QGridLayout(box)
        layout.setHorizontalSpacing(10)
        layout.setVerticalSpacing(10)

        self.input_entry = QLineEdit()
        self.input_entry.setPlaceholderText("Select an .xlsx input file")
        self.input_entry.setToolTip("Select the Excel input file containing taxa data.")

        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self.browse_input_file)

        self.output_entry = QLineEdit()
        self.output_entry.setPlaceholderText("Enter project name for results and config")
        self.output_entry.setToolTip("Defines the result directory and output file names.")

        output_hint = QLabel("Output folder and files will be created in the root directory.")
        output_hint.setWordWrap(True)
        output_hint.setObjectName("infoHint")

        layout.addWidget(QLabel("Input File"), 0, 0)
        layout.addWidget(self.input_entry, 0, 1)
        layout.addWidget(browse_btn, 0, 2)
        layout.addWidget(QLabel("Output Name"), 1, 0)
        layout.addWidget(self.output_entry, 1, 1, 1, 2)
        layout.addWidget(output_hint, 2, 0, 1, 3)
        return box

    def create_taxa_box(self) -> QGroupBox:
        box = self.section_box("Options")
        layout = QGridLayout(box)
        layout.setHorizontalSpacing(14)
        layout.setVerticalSpacing(10)

        self.taxa_col_entry = QLineEdit("1")
        self.taxa_col_entry.setToolTip("Column number containing taxa names.")
        self.taxa_col_entry.setFixedWidth(55)

        self.header_row_entry = QLineEdit("1")
        self.header_row_entry.setToolTip("Row number containing the table header.")
        self.header_row_entry.setFixedWidth(55)

        self.id_column_entry = QLineEdit("false")
        self.id_column_entry.setToolTip("Column containing IDs for taxa assignment.")
        self.id_column_entry.setFixedWidth(70)

        layout.addWidget(QLabel("Taxa Column"), 0, 0)
        layout.addWidget(self.taxa_col_entry, 0, 1)
        layout.addWidget(QLabel("Header Row"), 0, 2)
        layout.addWidget(self.header_row_entry, 0, 3)
        layout.addWidget(QLabel("ID Column"), 0, 4)
        layout.addWidget(self.id_column_entry, 0, 5)
        layout.setColumnStretch(6, 1)
        return box

    def create_options_box(self) -> QGroupBox:
        box = self.section_box("Options")
        layout = QGridLayout(box)
        layout.setHorizontalSpacing(14)
        layout.setVerticalSpacing(8)

        self.backup_config_check = QCheckBox("Backup Config")
        self.backup_input_check = QCheckBox("Backup Input")
        self.export_missing_check = QCheckBox("Export Missing")
        self.taxa_fallback_check = QCheckBox("Taxa Fallback")
        self.backup_blacklist_check = QCheckBox("Backup Blacklist")
        self.backup_skiplist_check = QCheckBox("Backup Skiplist")

        self.backup_config_check.setChecked(True)
        self.backup_input_check.setChecked(True)
        self.export_missing_check.setChecked(True)
        self.taxa_fallback_check.setChecked(False)
        self.backup_blacklist_check.setChecked(True)
        self.backup_skiplist_check.setChecked(True)

        layout.addWidget(self.backup_config_check, 0, 0)
        layout.addWidget(self.backup_input_check, 0, 1)
        layout.addWidget(self.export_missing_check, 0, 2)
        layout.addWidget(self.taxa_fallback_check, 1, 0)
        layout.addWidget(self.backup_blacklist_check, 1, 1)
        layout.addWidget(self.backup_skiplist_check, 1, 2)
        return box

    def create_database_box(self) -> QGroupBox:
        box = self.section_box("Databases")
        layout = QVBoxLayout(box)
        layout.setSpacing(10)

        algaebase_hint = QLabel("AlgaeBase DB requires an API key!")
        algaebase_hint.setWordWrap(True)
        algaebase_hint.setObjectName("warningHint")

        checks = QHBoxLayout()
        self.ncbi_check = QCheckBox("NCBI")
        self.pr2_check = QCheckBox("PR2")
        self.algaebase_check = QCheckBox("AlgaeBase")
        self.ncbi_check.setChecked(True)
        self.pr2_check.setChecked(True)
        self.algaebase_check.setChecked(True)

        checks.addWidget(self.ncbi_check)
        checks.addWidget(self.pr2_check)
        checks.addWidget(self.algaebase_check)
        checks.addStretch()

        layout.addWidget(algaebase_hint)
        layout.addSpacing(5)
        layout.addLayout(checks)
        return box

    def create_pr2_box(self) -> QGroupBox:
        box = self.section_box("PR2 Database")
        layout = QGridLayout(box)
        layout.setHorizontalSpacing(10)
        layout.setVerticalSpacing(10)

        pr2_hint = QLabel("Older database versions can also be used.")
        pr2_hint.setWordWrap(True)
        pr2_hint.setObjectName("infoHint")

        self.pr2_entry = QLineEdit(DEFAULT_PR2_PATH)
        self.pr2_entry.setCursorPosition(0)
        self.pr2_entry.setToolTip("Path to the local PR2 taxonomy database file.")

        pr2_browse_btn = QPushButton("Browse")
        pr2_browse_btn.clicked.connect(self.browse_pr2_database)

        layout.addWidget(pr2_hint, 0, 0, 1, 3)
        layout.addWidget(QLabel("PR2 Database"), 1, 0)
        layout.addWidget(self.pr2_entry, 1, 1)
        layout.addWidget(pr2_browse_btn, 1, 2)
        return box

    def create_ncbi_box(self) -> QGroupBox:
        box = self.section_box("NCBI API")
        layout = QGridLayout(box)
        layout.setHorizontalSpacing(10)
        layout.setVerticalSpacing(10)

        self.email_entry = QLineEdit(DEFAULT_NCBI_EMAIL)
        self.email_entry.setToolTip("Email address used for NCBI API requests.")

        self.api_key_entry = QLineEdit()
        self.api_key_entry.setPlaceholderText("Optional API key for higher NCBI limits.")
        self.api_key_entry.setToolTip("Optional, but recommended for higher NCBI request limits.")
        self.api_key_entry.setEchoMode(QLineEdit.EchoMode.Normal)

        self.hide_ncbi_check = QCheckBox("Hide")
        self.hide_ncbi_check.stateChanged.connect(
            lambda: self.toggle_password_visibility(self.api_key_entry, self.hide_ncbi_check)
        )

        layout.addWidget(QLabel("NCBI Email"), 0, 0)
        layout.addWidget(self.email_entry, 0, 1, 1, 2)
        layout.addWidget(QLabel("NCBI API Key"), 1, 0)
        layout.addWidget(self.api_key_entry, 1, 1)
        layout.addWidget(self.hide_ncbi_check, 1, 2)
        return box

    def create_algaebase_box(self) -> QGroupBox:
        box = self.section_box("AlgaeBase API")
        layout = QGridLayout(box)
        layout.setHorizontalSpacing(10)
        layout.setVerticalSpacing(10)

        self.algb_key_entry = QLineEdit()
        self.algb_key_entry.setPlaceholderText("Required AlgaeBase API key")
        self.algb_key_entry.setToolTip("Required when AlgaeBase is enabled.")

        self.algb_url_entry = QLineEdit(DEFAULT_ALGAEBASE_URL)
        self.algb_url_entry.setToolTip("AlgaeBase API endpoint used for species queries.")

        self.hide_algb_check = QCheckBox("Hide")
        self.hide_algb_check.stateChanged.connect(
            lambda: self.toggle_password_visibility(self.algb_key_entry, self.hide_algb_check)
        )

        layout.addWidget(QLabel("ALGB API Key"), 0, 0)
        layout.addWidget(self.algb_key_entry, 0, 1)
        layout.addWidget(self.hide_algb_check, 0, 2)
        layout.addWidget(QLabel("ALGB API URL"), 1, 0)
        layout.addWidget(self.algb_url_entry, 1, 1, 1, 2)
        return box

    def create_controls_box(self) -> QGroupBox:
        box = self.section_box("Controls")
        layout = QVBoxLayout(box)
        layout.setSpacing(12)

        control_hint = QLabel("Save the configuration before running the query!")
        control_hint.setWordWrap(True)
        control_hint.setObjectName("warningHint")

        self.process_status = QLabel("Status: Idle")
        self.process_status.setWordWrap(True)
        self.process_status.setObjectName("processIdle")

        save_btn = QPushButton("Save Config")
        start_btn = QPushButton("Start Query")
        load_btn = QPushButton("Load Config")
        reset_btn = QPushButton("Reset")

        save_btn.clicked.connect(self.save_config)
        start_btn.clicked.connect(self.start_query)
        load_btn.clicked.connect(self.load_config)
        reset_btn.clicked.connect(self.reset_fields)

        button_row = QHBoxLayout()
        button_row.setSpacing(10)
        button_row.addWidget(save_btn)
        button_row.addWidget(start_btn)
        button_row.addWidget(load_btn)
        button_row.addWidget(reset_btn)
        button_row.addStretch()

        layout.addWidget(control_hint)
        layout.addWidget(self.process_status)
        layout.addLayout(button_row)
        return box

    def section_box(self, title: str) -> QGroupBox:
        """Create a styled group box used for all GUI sections."""
        box = QGroupBox(title)
        box.setObjectName("sectionBox")
        return box

    def browse_input_file(self):
        """Open a file dialog and store the selected input Excel file path."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Excel file",
            str(Path.cwd()),
            "Excel file (*.xlsx)",
        )
        if file_path:
            self.input_entry.setText(file_path)

    def browse_pr2_database(self):
        """Open a file dialog and store the selected local PR2 database path."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select PR2 database file",
            str(Path.cwd()),
            "Excel file (*.xlsx)",
        )
        if file_path:
            self.pr2_entry.setText(file_path)

    def toggle_password_visibility(self, entry: QLineEdit, checkbox: QCheckBox):
        """Hide or show API key fields based on the matching checkbox."""
        entry.setEchoMode(
            QLineEdit.EchoMode.Password if checkbox.isChecked() else QLineEdit.EchoMode.Normal
        )

    def set_process_status(self, text: str, style_name: str):
        """Update process status text and refresh its Qt stylesheet state."""
        self.process_status.setText(text)
        self.process_status.setObjectName(style_name)
        self.process_status.style().unpolish(self.process_status)
        self.process_status.style().polish(self.process_status)

    def yaml_bool(self, value: bool) -> str:
        return "true" if value else "false"

    def yaml_quote(self, value: str) -> str:
        escaped_value = value.replace("\\", "/").replace('"', '\\"')
        return f'"{escaped_value}"'

    def yaml_number_or_false(self, value: str) -> str:
        cleaned_value = value.strip()
        if cleaned_value.lower() == "false":
            return "false"
        if cleaned_value.isdigit():
            return cleaned_value
        return self.yaml_quote(cleaned_value)

    def config_file_path(self, output_name: str) -> Path:
        """Return a safe config file path derived from the output/project name."""
        safe_name = "".join(
            char if char.isalnum() or char in ("_", "-", ".") else "_"
            for char in output_name.strip()
        ).strip("._")

        config_path = Path.cwd() / (safe_name or "algaetax_config")
        if config_path.suffix.lower() not in {".yaml", ".yml"}:
            config_path = config_path.with_suffix(".yaml")
        return config_path

    def build_config_text(self) -> str:
        """Create the YAML configuration text from all current GUI values."""
        input_data = self.input_entry.text().strip()
        output_name = self.output_entry.text().strip()
        taxa_column_number = self.yaml_number_or_false(self.taxa_col_entry.text())
        header_row = self.yaml_number_or_false(self.header_row_entry.text())
        id_column_number = self.yaml_number_or_false(self.id_column_entry.text())

        backup_config = self.yaml_bool(self.backup_config_check.isChecked())
        backup_input = self.yaml_bool(self.backup_input_check.isChecked())
        export_missing = self.yaml_bool(self.export_missing_check.isChecked())
        synonym_fallback = self.yaml_bool(self.taxa_fallback_check.isChecked())
        backup_blacklist = self.yaml_bool(self.backup_blacklist_check.isChecked())
        backup_skiplist = self.yaml_bool(self.backup_skiplist_check.isChecked())

        use_ncbi = self.yaml_bool(self.ncbi_check.isChecked())
        use_pr2 = self.yaml_bool(self.pr2_check.isChecked())
        use_algaebase = self.yaml_bool(self.algaebase_check.isChecked())

        pr2_database = self.pr2_entry.text().strip()
        ncbi_api_key = self.api_key_entry.text().strip()
        ncbi_email = self.email_entry.text().strip()
        algaebase_api_key = self.algb_key_entry.text().strip()
        algaebase_api_url = self.algb_url_entry.text().strip()

        return f"""# algaetax; configuration file: config.yaml

# General settings
general:
  input_data: {self.yaml_quote(input_data)} # Path to input Excel file (e.g. input_data/example_file.xlsx)
  output_dir: {self.yaml_quote(output_name)} # Folder to save results (default: results_algaetax/)
  taxa_column_number: {taxa_column_number} # 1-based index: 1 = column A, 27 = column AA
  header_row: {header_row} # Excel header row (1 = first row); set to false if the file has no header
  backup_config: {backup_config} # Copy the config.yaml file to the output directory for reference (default=true)
  backup_input: {backup_input} # Backup of the input Excel file before processing (true/false)
  export_not_found_taxa_list: {export_missing} # Export list of taxa where ALL used databases returned "Not found" (file: taxa_not_found.csv)
  synonym_fallback: {synonym_fallback} # Enable synonym fallback via AlgaeBase for NCBI/PR2 if not found

  # Optional: defines which column in the input file contains unique sample or sequence IDs.
  # These IDs will be added to the output for easier tracking of each taxon entry.
  id_column_number: {id_column_number}  # Use 1-based index (e.g., 1 = column A); set to false to disable

# Taxa filter settings
filter:
  blacklist_file: "blacklist.txt" # Path to text file with terms to ignore in species names
  backup_blacklist: {backup_blacklist} # Create a backup of the blacklist in {{output_dir}}/backups (default=true)
  skiplist_file: "skiplist.txt" # Terms that cause taxa to be skipped from all DB queries but kept in output
  backup_skiplist: {backup_skiplist} # Save a copy of skiplist.txt into the {{output_dir}}/backups folder

# Databases to use
database:
  NCBI: {use_ncbi} # Use NCBI Database
  PR2: {use_pr2} # Use local PR2 Database
  ALGB: {use_algaebase} # Use AlgaeBase Database (requires API key!)

# Paths to local database files
database_path:
  db_pr2: {self.yaml_quote(pr2_database)} # Path to PR2 file (e.g. database/pr2_v5.xlsx)

# NCBI API settings
db_ncbi:
  api_key: {self.yaml_quote(ncbi_api_key)} # Optional key to speed up requests (e.g. abc123xyz)
  ncbi_email: {self.yaml_quote(ncbi_email)} # Email for NCBI queries (default: anonymous@ncbi.com, e.g. name@example.com)

# AlgaeBase API settings
db_algaebase:
  api_key: {self.yaml_quote(algaebase_api_key)} # API key for AlgaeBase (e.g. abc123xyz)
  api_url: {self.yaml_quote(algaebase_api_url)} # API URL (see docs: algaebase.org/api/)
"""

    def save_config(self):
        """Validate required values and write the YAML configuration file."""
        output_name = self.output_entry.text().strip()

        if self.algaebase_check.isChecked() and not self.algb_key_entry.text().strip():
            QMessageBox.warning(
                self,
                "Missing AlgaeBase API Key",
                "AlgaeBase is enabled but no API key was provided. Add an API key or disable AlgaeBase.",
            )
            return

        if not output_name:
            QMessageBox.warning(
                self,
                "Missing Output Name",
                "Please enter an output name before saving. Add an output name or cancel the configuration save.",
            )
            return

        config_path = self.config_file_path(output_name)
        if config_path.exists() and not self.confirm_overwrite_config(config_path):
            return

        try:
            config_path.write_text(self.build_config_text(), encoding="utf-8")
        except OSError as error:
            QMessageBox.critical(self, "Save Config", f"Could not save configuration file:\n{error}")
            return

        QMessageBox.information(self, "Save Config", f"Configuration saved as:\n{config_path}")

    def confirm_overwrite_config(self, config_path: Path) -> bool:
        """Ask the user before overwriting an existing config file."""
        reply = QMessageBox.question(
            self,
            "Overwrite Config",
            f"The configuration file already exists:\n\n{config_path.name}\n\nDo you want to overwrite it?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        return reply == QMessageBox.StandardButton.Yes

    def start_query(self):
        """Start algaetax in an external terminal using the saved config file."""
        output_name = self.output_entry.text().strip()
        if not output_name:
            QMessageBox.warning(self, "Missing Output Name", "Please enter an output name before starting the query.")
            return

        config_path = self.config_file_path(output_name)
        output_dir = Path.cwd() / output_name

        if output_dir.exists() and not self.confirm_existing_output_dir(output_dir):
            return

        if not config_path.exists():
            QMessageBox.warning(self, "Missing Config File", "Config file not found. Please save the configuration first.")
            return

        if not (Path.cwd() / "algaetax.py").exists():
            QMessageBox.critical(self, "Missing Script", "algaetax.py was not found in the root directory.")
            return

        self.process_status_file = self.create_process_status_path(output_name)
        self.process_status_file.unlink(missing_ok=True)

        command = self.build_query_command(config_path)
        self.set_process_status("Status: Query in progress...", "processRunning")
        self.process_timer.start()

        self.process = QProcess(self)
        self.process.setWorkingDirectory(str(Path.cwd()))
        self.process.start("gnome-terminal", ["--", "bash", "-c", command])

        if not self.process.waitForStarted(3000):
            self.process_timer.stop()
            self.set_process_status("Status: Unable to launch terminal.", "processFailed")
            QMessageBox.critical(self, "Start Query", "Could not open terminal.")

    def confirm_existing_output_dir(self, output_dir: Path) -> bool:
        """Ask the user before running a query that may overwrite output files."""
        reply = QMessageBox.question(
            self,
            "Existing Output Directory",
            (
                f"The output directory already exists:\n\n{output_dir.name}\n\n"
                "Existing files may be overwritten.\nDo you want to continue?"
            ),
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        return reply == QMessageBox.StandardButton.Yes

    def create_process_status_path(self, output_name: str) -> Path:
        """Create a hidden status-file path used to track the external process."""
        safe_status_name = "".join(
            char if char.isalnum() or char in ("_", "-", ".") else "_"
            for char in output_name
        ).strip("._") or "algaetax"
        return Path.cwd() / f".{safe_status_name}_process_status.txt"

    def build_query_command(self, config_path: Path) -> str:
        """Build the shell command executed in the external terminal."""
        return (
            f'python3 algaetax.py --configfile "{config_path.name}"; '
            'exit_code=$?; '
            f'echo "$exit_code" > "{self.process_status_file.name}"; '
            'echo ""; '
            'if [ "$exit_code" -eq 0 ]; then '
            'echo "Analysis finished successfully."; '
            'else '
            'echo "Analysis failed."; '
            'fi; '
            'read -p "Press Enter to close..."; '
            'exit "$exit_code"'
        )

    def check_process_status(self):
        """Poll the status file to detect when the external query has finished."""
        if not self.process_status_file or not self.process_status_file.exists():
            return

        exit_code = self.process_status_file.read_text(encoding="utf-8").strip()
        self.process_status_file.unlink(missing_ok=True)
        self.process_timer.stop()

        if exit_code == "0":
            self.set_process_status("Status: Query completed successfully.", "processSuccess")
        else:
            self.set_process_status("Status: Query failed.", "processFailed")

    def reset_fields(self):
        """Restore all GUI fields to their default values."""
        self.input_entry.clear()
        self.output_entry.clear()
        self.taxa_col_entry.setText("1")
        self.header_row_entry.setText("1")
        self.id_column_entry.setText("false")
        self.backup_config_check.setChecked(True)
        self.backup_input_check.setChecked(True)
        self.export_missing_check.setChecked(True)
        self.taxa_fallback_check.setChecked(False)
        self.backup_blacklist_check.setChecked(True)
        self.backup_skiplist_check.setChecked(True)
        self.ncbi_check.setChecked(True)
        self.pr2_check.setChecked(True)
        self.algaebase_check.setChecked(True)
        self.pr2_entry.setText(DEFAULT_PR2_PATH)
        self.pr2_entry.setCursorPosition(0)
        self.email_entry.setText(DEFAULT_NCBI_EMAIL)
        self.api_key_entry.clear()
        self.algb_key_entry.clear()
        self.algb_url_entry.setText(DEFAULT_ALGAEBASE_URL)
        self.set_process_status("Status: Idle", "processIdle")

    def strip_yaml_comment(self, value: str) -> str:
        """Remove inline YAML comments while preserving quoted text."""
        in_quotes = False
        quote_char = ""

        for index, char in enumerate(value):
            if char in {"'", '"'}:
                if not in_quotes:
                    in_quotes = True
                    quote_char = char
                elif quote_char == char:
                    in_quotes = False

            if char == "#" and not in_quotes:
                return value[:index].strip()

        return value.strip()

    def yaml_unquote(self, value: str) -> str:
        """Remove surrounding YAML quotes before loading values into the GUI."""
        cleaned_value = self.strip_yaml_comment(value).strip()
        if (
            len(cleaned_value) >= 2
            and cleaned_value[0] == cleaned_value[-1]
            and cleaned_value[0] in {"'", '"'}
        ):
            cleaned_value = cleaned_value[1:-1]
        return cleaned_value.replace('\\"', '"')

    def yaml_to_bool(self, value: str) -> bool:
        return self.yaml_unquote(value).lower() == "true"

    def parse_config_file(self, file_path: Path) -> dict:
        """Parse the simple YAML structure generated by this GUI."""
        config = {}
        current_section = None

        for line in file_path.read_text(encoding="utf-8").splitlines():
            stripped_line = line.strip()

            if not stripped_line or stripped_line.startswith("#"):
                continue

            if not line.startswith(" ") and stripped_line.endswith(":"):
                current_section = stripped_line[:-1]
                config[current_section] = {}
                continue

            if current_section and ":" in stripped_line:
                key, value = stripped_line.split(":", 1)
                config[current_section][key.strip()] = self.strip_yaml_comment(value)

        return config

    def load_config(self):
        """Load an existing YAML configuration file into the GUI fields."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select YAML file",
            str(Path.cwd()),
            "YAML file (*.yaml *.yml)",
        )
        if not file_path:
            return

        config_path = Path(file_path)
        try:
            config = self.parse_config_file(config_path)
        except OSError as error:
            QMessageBox.critical(self, "Load Config", f"Could not load configuration file:\n{error}")
            return

        general = config.get("general", {})
        filter_settings = config.get("filter", {})
        database = config.get("database", {})
        database_path = config.get("database_path", {})
        db_ncbi = config.get("db_ncbi", {})
        db_algaebase = config.get("db_algaebase", {})

        self.input_entry.setText(self.yaml_unquote(general.get("input_data", "")))
        self.output_entry.setText(self.yaml_unquote(general.get("output_dir", "")))
        self.taxa_col_entry.setText(self.yaml_unquote(general.get("taxa_column_number", "1")))
        self.header_row_entry.setText(self.yaml_unquote(general.get("header_row", "1")))
        self.id_column_entry.setText(self.yaml_unquote(general.get("id_column_number", "false")))

        self.backup_config_check.setChecked(self.yaml_to_bool(general.get("backup_config", "true")))
        self.backup_input_check.setChecked(self.yaml_to_bool(general.get("backup_input", "true")))
        self.export_missing_check.setChecked(self.yaml_to_bool(general.get("export_not_found_taxa_list", "true")))
        self.taxa_fallback_check.setChecked(self.yaml_to_bool(general.get("synonym_fallback", "false")))
        self.backup_blacklist_check.setChecked(self.yaml_to_bool(filter_settings.get("backup_blacklist", "true")))
        self.backup_skiplist_check.setChecked(self.yaml_to_bool(filter_settings.get("backup_skiplist", "true")))
        self.ncbi_check.setChecked(self.yaml_to_bool(database.get("NCBI", "true")))
        self.pr2_check.setChecked(self.yaml_to_bool(database.get("PR2", "true")))
        self.algaebase_check.setChecked(self.yaml_to_bool(database.get("ALGB", "true")))

        self.pr2_entry.setText(self.yaml_unquote(database_path.get("db_pr2", DEFAULT_PR2_PATH)))
        self.pr2_entry.setCursorPosition(0)
        self.api_key_entry.setText(self.yaml_unquote(db_ncbi.get("api_key", "")))
        self.email_entry.setText(self.yaml_unquote(db_ncbi.get("ncbi_email", DEFAULT_NCBI_EMAIL)))
        self.algb_key_entry.setText(self.yaml_unquote(db_algaebase.get("api_key", "")))
        self.algb_url_entry.setText(self.yaml_unquote(db_algaebase.get("api_url", DEFAULT_ALGAEBASE_URL)))

        QMessageBox.information(self, "Load Config", f"Configuration loaded from:\n{config_path}")

    def apply_style(self):
        """Apply the dark theme and custom widget styling."""
        self.setStyleSheet(
            """
            QWidget {
                background-color: #202323;
                color: #e2e5e4;
                font-family: Arial;
                font-size: 15px;
            }

            AlgaeTaxGUI {
                background-color: #1C1C1C;
                border: 2px solid #43514f;
            }

            QLabel {
                background-color: transparent;
                color: #e2e5e4;
            }

            QLabel#title {
                font-size: 25px;
                color: #f1f3f2;
            }

            QGroupBox#sectionBox {
                background-color: #1F1F1F;
                border: 2px solid #5F6E5D;
                border-radius: 6px;
                margin-top: 14px;
                padding: 16px;
                color: #78aaa3;
            }

            QGroupBox#sectionBox::title {
                subcontrol-origin: margin;
                left: 14px;
                padding: 6px 8px;
                color: #e5e5e5;
                font-size: 18px;
            }

            QLineEdit {
                background-color: #181b1b;
                color: #f0f2f1;
                border: 2px solid #4a5553;
                border-radius: 4px;
                padding: 7px 9px;
                min-height: 24px;
                selection-background-color: #59706C;
                selection-color: #ffffff;
            }

            QLineEdit:focus {
                border: 2px solid #78aaa3;
                background-color: #1d2020;
            }

            QPushButton {
                background-color: #52756F;
                color: #f4f6f5;
                border: 2px solid #52756F;
                border-radius: 4px;
                padding: 10px 17px;
            }

            QPushButton:hover {
                background-color: #657f7a;
                border: 2px solid #86a7a1;
            }

            QPushButton:pressed {
                background-color: #455956;
                border: 2px solid #78aaa3;
            }

            QCheckBox {
                spacing: 10px;
                background-color: transparent;
                color: #e2e5e4;
            }

            QCheckBox::indicator {
                width: 22px;
                height: 22px;
            }

            QCheckBox:disabled {
                color: #777f7d;
            }

            QMessageBox {
                background-color: #282c2c;
                color: #e2e5e4;
            }

            QLabel#infoHint {
                color: #f1f5f3;
                background-color: #313938;
                border-radius: 4px;
                padding: 7px 9px;
            }

            QLabel#warningHint {
                color: #f1f5f3;
                background-color: #3E5943;
                border-radius: 4px;
                padding: 7px 9px;
            }

            QLabel#processIdle,
            QLabel#processRunning,
            QLabel#processSuccess,
            QLabel#processFailed {
                color: #ffffff;
                background-color: #181b1b;
                border: 1px solid #4a5553;
                border-radius: 4px;
                padding: 7px 9px;
            }
            """
        )


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle(DarkCheckboxStyle())
    window = AlgaeTaxGUI()
    window.show()
    sys.exit(app.exec())
