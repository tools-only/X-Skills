# claude-openocd-spi-dump

**Research Date**: 2026-02-20
**Source URL**: <https://github.com/lukejenkins/claude-openocd-spi-dump>
**GitHub Repository**: <https://github.com/lukejenkins/claude-openocd-spi-dump>
**Version at Research**: 1.0.0
**License**: No license file (public repository)

---

## Overview

A Claude Code plugin that enables AI-assisted dumping of SPI flash/EEPROM memory through a microcontroller's debug interface (SWD/JTAG) via OpenOCD, without requiring dedicated SPI programming hardware. The plugin provides comprehensive MCU register knowledge, code templates, and guided interactive workflows for firmware extraction from embedded systems.

---

## Problem Addressed

| Problem | Solution |
|---------|----------|
| Extracting firmware from SPI flash requires expensive dedicated programming hardware (CH341A, Bus Pirate, etc.) | Leverages existing debug interfaces (SWD/JTAG) that developers already have for reverse engineering work |
| SPI flash chips on production boards often lack accessible test points or are in BGA packages | Reads flash indirectly through the microcontroller's SPI peripheral using RAM-resident code |
| Writing low-level embedded code requires extensive MCU datasheet research and register-level programming knowledge | Provides ready-to-use register maps for 6 MCU families (SAM4S, SAM3X, STM32F1/F4, nRF52, LPC1768) and customizable templates |
| Debugging SPI communication issues and hardware setup requires deep embedded systems expertise | Offers guided troubleshooting with heartbeat monitoring, JEDEC ID verification, and comprehensive error diagnostics |
| Embedded developers need to context-switch between datasheets, code examples, and debugging tools | Integrates all knowledge into Claude's context enabling conversational guidance through the entire process |

---

## Key Statistics

| Metric | Value | Date Gathered |
|--------|-------|---------------|
| GitHub Stars | 3 | 2026-02-20 |
| Forks | 0 | 2026-02-20 |
| Open Issues | 0 | 2026-02-20 |
| Primary Language | Shell | 2026-02-20 |
| Created | 2025-12-25 | 2026-02-20 |
| Last Push | 2025-12-28 | 2026-02-20 |
| Topics | claude-code, claude-code-plugin, embedded, firmware-extraction, openocd, reverse-engineering, spi-flash | 2026-02-20 |

---

## Key Features

### Domain-Specific Knowledge Integration

- **MCU Register Maps**: Complete register definitions for SPI, GPIO, watchdog, and vector table configuration for 6 microcontroller families
- **SPI Flash Commands**: JEDEC command reference (READ 0x03, RDID 0x9F, RDSR 0x05) with manufacturer ID decoding
- **Memory Layout Planning**: SRAM allocation strategies for vector tables, code, buffers, stack, and communication areas
- **Critical Implementation Details**: Documented gotchas like VTOR initialization, SAM4S PCS field requirements, GPIO CS reclamation

### RAM-Resident Code Architecture

- **Minimal Footprint**: ~500 byte programs that execute entirely from MCU SRAM without flash dependency
- **Vector Table Handling**: Proper Cortex-M vector table setup with VTOR configuration for exception handling
- **Host-MCU Communication Protocol**: Memory-mapped command/status interface for coordinating operations between OpenOCD and running code
- **Heartbeat Monitoring**: Real-time execution verification to detect hung states and polling loops

### Guided Interactive Workflow

- **`/spi-dump` Command**: Multi-phase workflow that gathers hardware information, generates customized code, and walks through testing
- **Adaptive Code Generation**: Customizes C source, linker scripts, and OpenOCD TCL based on user's specific MCU and SPI configuration
- **Progressive Testing**: Starts with JEDEC ID verification before full flash dumps to catch configuration issues early
- **Troubleshooting Integration**: Context-aware error diagnosis based on symptoms (HardFault, SPI hangs, CS issues)

### Automation Scripts

- **MCU-Agnostic Shell Scripts**: `init_dumper.sh`, `dump.sh`, `verify_dump.sh` work across all supported MCU families
- **Proper Initialization**: Reads SP/PC from vector table automatically to ensure code starts from reset vector
- **Parameterized Memory Addresses**: Environment variables (SRAM_BASE, COMM_BASE) for different MCU families
- **Dump Verification**: Integrity checking with stuck data line detection

### Extensibility

- **Clear Adaptation Guide**: Step-by-step instructions for adding new MCU families
- **Template-Based Design**: Reference implementations demonstrate patterns applicable to any Cortex-M device
- **Portable Algorithm**: Core SPI dumping logic remains constant across MCUs, only register addresses change

---

## Technical Architecture

### Three-Layer Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                      Claude Code Layer                      │
│  - Skill: Comprehensive MCU/SPI knowledge and workflows     │
│  - Command: Interactive guided setup and code generation   │
└────────────────────┬────────────────────────────────────────┘
                     │
┌────────────────────▼────────────────────────────────────────┐
│                     Host PC Layer                           │
│  - OpenOCD: Debug interface, memory read/write, TCL scripts│
│  - Build Tools: arm-none-eabi-gcc for RAM-resident binary   │
│  - Shell Scripts: Automation for load/dump/verify cycle    │
└────────────────────┬────────────────────────────────────────┘
                     │ SWD/JTAG
┌────────────────────▼────────────────────────────────────────┐
│                    Target MCU Layer                         │
│  - SRAM: Executes dumper code loaded via debug interface   │
│  - SPI Peripheral: Controlled via register writes          │
│  - GPIO: CS pin toggling for flash chip select             │
└────────────────────┬────────────────────────────────────────┘
                     │ SPI Bus
┌────────────────────▼────────────────────────────────────────┐
│                   SPI Flash/EEPROM                          │
│  - Responds to standard JEDEC commands                      │
│  - Returns data through MCU SPI peripheral                  │
└─────────────────────────────────────────────────────────────┘
```

### Memory Layout Strategy

```
SRAM (example: 64KB at 0x20000000)

0x20000000  ┌─────────────────┐
            │ Vector Table    │  64 bytes (16 Cortex-M exception vectors)
            │ SP: 0x2000FF00  │  Initial stack pointer
            │ PC: Reset vector│  Entry point address
0x20000040  ├─────────────────┤
            │ Dumper Code     │  ~400-600 bytes
            │ - spi_init()    │  Initialize SPI peripheral
            │ - spi_transfer()│  Byte-level SPI transaction
            │ - command_loop()│  Poll comm area, execute commands
0x2000E000  ├─────────────────┤
            │ Read Buffer     │  4KB (tunable based on SRAM size)
            │                 │  Flash data staged here before host read
0x2000FE00  ├─────────────────┤
            │ Stack Space     │  256 bytes (grows downward)
0x2000FF00  ├─────────────────┤
            │ Communication   │  256 bytes
            │ Area            │  STATUS, FLASH_ADDR, SIZE, DEST,
            │                 │  JEDEC_ID, ERROR, HEARTBEAT
0x20010000  └─────────────────┘
```

### Communication Protocol

Host (OpenOCD) and target (MCU code) communicate via memory-mapped registers:

| Offset | Register | Direction | Purpose |
|--------|----------|-----------|---------|
| +0x00 | STATUS | Bidirectional | Command from host (0x10=Read, 0x20=JEDEC, 0xFF=Exit), status to host (0x00=Idle, 0x01=Busy, 0x02=Done, 0xDEADxxxx=Error) |
| +0x04 | FLASH_ADDR | Host → MCU | 24-bit SPI flash address to read |
| +0x08 | SIZE | Host → MCU | Number of bytes to read |
| +0x0C | DEST | Host → MCU | Destination buffer address in SRAM |
| +0x10 | JEDEC_ID | MCU → Host | Result of 0x9F JEDEC ID read (3 bytes: Manufacturer, Type, Capacity) |
| +0x14 | ERROR | MCU → Host | Error code details when STATUS indicates error |
| +0x18 | HEARTBEAT | MCU → Host | Increments in main loop to prove execution |

### Execution Flow

1. **Build Phase**: arm-none-eabi-gcc compiles C source with custom linker script placing all code/data in SRAM
2. **Load Phase**: OpenOCD writes binary to target SRAM via SWD/JTAG
3. **Init Phase**: OpenOCD reads SP and PC from vector table, sets registers, resumes execution
4. **Verify Phase**: Host reads JEDEC ID via 0x20 command to validate SPI communication
5. **Dump Phase**: Host issues 0x10 commands with address/size, reads buffer after MCU sets STATUS=Done
6. **Repeat**: Chunk through entire flash address space until complete

---

## Installation & Usage

### Prerequisites

**Required Software:**

```bash
# macOS
brew install openocd arm-none-eabi-gcc

# Linux (Debian/Ubuntu)
apt install openocd gcc-arm-none-eabi
```

**Required Hardware:**
- SWD/JTAG debug probe (ST-Link, J-Link, CMSIS-DAP)
- Working OpenOCD connection to target MCU
- Target MCU with SPI flash/EEPROM connected to SPI peripheral

### Installation

```bash
# Use as plugin directory
claude --plugin-dir /path/to/claude-openocd-spi-dump

# Or copy to project's .claude-plugin directory
cp -r claude-openocd-spi-dump /path/to/project/.claude-plugin/
```

### Quick Start with Guided Workflow

```bash
# In Claude Code session with plugin loaded
/spi-dump
```

This launches an interactive workflow that:
1. Gathers MCU and SPI flash information
2. Generates customized C code, linker script, and OpenOCD TCL
3. Provides build instructions
4. Guides through JEDEC ID verification
5. Orchestrates full flash dump

### Manual Workflow with Shell Scripts

```bash
# 1. Customize spi_dump.c for your MCU (update register addresses)
# 2. Build the dumper
arm-none-eabi-gcc -mcpu=cortex-m4 -mthumb -Os -ffreestanding \
    -nostdlib -T spi_dump.ld -o spi_dump.elf spi_dump.c
arm-none-eabi-objcopy -O binary spi_dump.elf spi_dump.bin

# 3. Initialize (automatically reads vector table)
./init_dumper.sh spi_dump.bin

# 4. Dump flash (4MB example)
./dump.sh firmware.bin 0x400000

# 5. Verify integrity
./verify_dump.sh firmware.bin 0x400000
```

**For LPC1768 (non-standard SRAM at 0x10000000):**

```bash
SRAM_BASE=0x10000000 ./init_dumper.sh spi_dump.bin
SRAM_BASE=0x10000000 ./dump.sh firmware.bin 0x100000
```

### Supported MCU Families

| Family | Examples | SRAM Base | Status |
|--------|----------|-----------|--------|
| SAM4S | ATSAM4S2A, ATSAM4S4A | 0x20000000 | Full register maps provided |
| SAM3X | ATSAM3X8E (Arduino Due) | 0x20000000 | Full register maps provided |
| STM32F1 | STM32F103 (Blue Pill) | 0x20000000 | Full register maps provided |
| STM32F4 | STM32F407, STM32F411 | 0x20000000 | Full register maps provided |
| nRF52 | nRF52832, nRF52840 | 0x20000000 | Full register maps provided |
| LPC1768 | LPC1768 (mbed) | 0x10000000 | Full register maps provided |

Adding new MCU families requires only finding SPI/GPIO register addresses in the datasheet—the core algorithm is portable.

---

## Relevance to Claude Code Development

### Applications

**AI-Assisted Hardware Reverse Engineering**: This plugin demonstrates how Claude Code can bridge software and hardware domains. By encoding deep embedded systems knowledge (register maps, peripheral initialization sequences, debugging strategies) into Claude's context, developers can perform complex hardware tasks through conversational interfaces without constantly referencing datasheets.

**Domain-Specific Tooling Pattern**: Shows effective strategy for bringing specialized technical domains into Claude Code ecosystem. The three-layer architecture (knowledge skill + guided command + automation scripts) is reusable for other hardware/firmware workflows like JTAG boundary scan, I2C EEPROM programming, or MCU peripheral configuration.

**Progressive Disclosure Workflow**: The `/spi-dump` command exemplifies how to structure complex multi-step technical processes as conversational workflows. Phases gather information, validate prerequisites, generate artifacts, and guide testing—a pattern applicable to any domain with setup/build/test cycles.

**Troubleshooting Knowledge Capture**: The comprehensive troubleshooting guide (symptom → cause → fix table) demonstrates how to encode expert debugging knowledge for AI retrieval. When users report "SPI hangs polling", Claude can immediately suggest checking SPI_MR PCS field based on documented patterns.

### Patterns Worth Adopting

**MCU Register Map as Structured Data**: Storing hardware register definitions in markdown tables enables both human readability and AI-assisted code generation. Claude can look up register addresses, generate initialization sequences, and explain hardware behavior without users consulting datasheets.

**Template + Customization Pipeline**: Providing complete working examples (C source, linker scripts, TCL) that get customized based on user context is more effective than explaining from scratch. Users get production-quality code adapted to their hardware.

**Heartbeat Monitoring for Embedded Debugging**: The heartbeat register pattern (incrementing counter proving execution) is a simple but powerful technique for debugging bare-metal code. Generalizable to any RAM-resident or early-boot code where traditional debugging isn't available.

**Memory-Mapped Communication Protocol**: Using fixed SRAM addresses for host-target communication is elegant and debugger-agnostic. This pattern works with any debug interface that can read/write memory (JTAG, SWD, BDM), not just OpenOCD.

**Verification-First Testing Strategy**: Requiring JEDEC ID success before full dump prevents wasting time on misconfigured hardware. Apply this "small test first" pattern to any hardware workflow (I2C probe before bulk read, GPIO toggle before complex protocol).

### Integration Opportunities

**Claude Code Skills for Other Embedded Tools**: This plugin validates the feasibility of similar tools for:
- JTAG boundary scan for board bring-up and testing
- I2C/SPI EEPROM programming and configuration
- MCU peripheral initialization (UART, ADC, PWM)
- Firmware patching and hot-loading for development

**Hardware Knowledge Base**: Could evolve into broader embedded systems knowledge plugin covering:
- Common MCU peripheral register patterns across vendors
- Debug interface protocols (SWD, JTAG, cJTAG)
- Bus protocol specifications (SPI, I2C, UART timing)
- Firmware binary format parsing (ELF, HEX, BIN)

**Integration with Ghidra/Binary Analysis**: Dumped firmware could flow directly into disassembly/reverse engineering workflows. Claude could orchestrate: dump → load in Ghidra → analyze → annotate → report findings.

**CI/CD for Embedded Systems**: The automation scripts (init, dump, verify) form foundation for firmware regression testing. Could expand to: build firmware → flash via debug → dump via SPI → compare against golden image.

**Educational Tool for Embedded Development**: The comprehensive documentation and guided workflows make this valuable for teaching:
- Cortex-M vector table and startup code
- SPI protocol implementation at register level
- OpenOCD TCL scripting
- RAM-resident code techniques

---

## References

- [GitHub Repository](https://github.com/lukejenkins/claude-openocd-spi-dump) (accessed 2026-02-20)
- [OpenOCD Documentation](http://openocd.org/documentation/) (accessed 2026-02-20)
- [ARM Cortex-M Programming Guide](https://developer.arm.com/documentation/dui0552/latest) (accessed 2026-02-20)
- Repository README.md (accessed 2026-02-20)
- skills/spi-flash-dump/SKILL.md (accessed 2026-02-20)
- commands/spi-dump.md (accessed 2026-02-20)

---

## Freshness Tracking

| Field | Value |
|-------|-------|
| Last Verified | 2026-02-20 |
| Version at Verification | 1.0.0 |
| Next Review Recommended | 2026-05-20 |

**Review Triggers:**
- New MCU family support added
- OpenOCD version compatibility changes
- Additional flash chip protocols supported
- Automation script enhancements
- GitHub star count increases significantly (>50 stars indicates growing adoption)
