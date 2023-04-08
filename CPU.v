// this is the code for CSC3050 project 4
// the code is written by Changyi Yang, 120090429
`timescale 1ns/1ns

module CPU();


    // generate the clock signal
    reg clock;
    wire on;
    parameter period = 50;
    assign on = 1;
    initial begin
    clock = 0;
    stall_times = 2'b00;
    #10000;
    forever begin
        #(period/2) clock = ~ clock;
    end
    // #(period/2) clock = ~ clock;
    // #(period/2) clock = ~ clock;
    // #(period/2) clock = ~ clock;
    // #(period/2) clock = ~ clock;
    // #(period/2) clock = ~ clock;
    // #(period/2) clock = ~ clock;
    // #(period/2) clock = ~ clock;
    // #(period/2) clock = ~ clock;
    // #(period/2) clock = ~ clock;
    // #(period/2) clock = ~ clock;
    // #(period/2) clock = ~ clock;
    // #(period/2) clock = ~ clock;


    // #100000000
    // $finish; 

    end 

    wire [31:0] instruction;

    // // the memory used to store instuction
    // reg [31:0] RAM [0:512-1];

    // // the registers
    // reg [31:0] Registers [31:0];

    // // the pipeline registers

    // // first part, some data needed to be passed


    reg [15:0] ALUControl_Reg;
    reg [32*4-1:0] SignImme_Reg;
    reg [3:0] ALUSrc_Reg;
    reg [3:0] MemtoReg_Reg;  
    reg [3:0] MemWrite_Reg;
    reg [2*4-1:0] Branch_Reg;
    reg [3:0] RegWrite_Reg;
    reg [2*4-1:0] RegDst_Reg;
    reg [3:0] shift_Reg;
    reg [32*4-1:0] PCF_Reg;
    initial PCF_Reg[32*4-1:32*3] = 32'h00000000;

    reg [32*4-1:0] rs_data_Reg; 
    reg [32*4-1:0] rt_data_Reg; 
    reg [5*4-1:0] rd_Reg;
    reg [5*4-1:0] sa_Reg;
    reg [32*4-1:0] MemOut_Reg;
    reg [32*4-1:0] PCPlus4_Reg;
    reg [32*4-1:0] ALUOut_Reg;
    reg [5*4-1:0] rt_Reg;

    reg [3:0] on_Reg;

    reg [3:0] jal_Reg;

    reg [1:0] stall_times;

    reg [2*4-1:0] Mem_WB_hazard_Reg;


// some wire to take the fresh-out output

    wire [31:0] PCF;
    wire jump;
    wire jump_src;
    wire [31:0] jump_address;
    wire [31:0] PCPlus4;


InstructionRAM instructionram(
    .clock(clock),
    .jump(jump),
    .jump_address(jump_address),
    .jump_src(jump_src),
    .rs_data(rs_data),

    .instruction(instruction),
    .PCPlus4(PCPlus4),

    .branch(branch),
    .branch_address(branch_address),

    .stall_times(lw_hazard),

    .jr_hazard(jr_hazard),
    .ALUOut(ALUOut),
    .ALUOut_2(ALUOut_Reg[32*2-1:32]),
    .ALUOut_3(ALUOut_Reg[32*3-1:32*2])
);



wire [31:0] rs_data;
wire [31:0] rt_data;
wire [4:0] rd;
wire [4:0] rt;
wire [3:0] ALUControl;
wire [4:0] sa;
wire shift;
wire [31:0] SignImme;
wire ALUSrc;
wire MemWrite;
wire MemtoReg;
wire RegWrite;
wire [1:0] RegDst;
wire on_sign;
wire slt_sign;
wire jal;
wire branch;
wire [31:0] branch_address;

wire [1:0] hazard_1;
wire [1:0] hazard_2;
wire [1:0] hazard_3;
wire [1:0] Mem_WB_hazard;
wire [1:0] lw_hazard;
wire [1:0] jr_hazard;


Register Reg(
    // the parameters for IF stage
    .clock(clock),
    .instruction_in(instruction),

    .rs_data(rs_data),
    .rt_data(rt_data),
    .rd_out(rd),
    .rt_out(rt),
    .ALUControl(ALUControl),
    .sa_out(sa),
    .shift(shift),
    .SignImme(SignImme),
    .ALUSrc(ALUSrc),
    .MemWrite(MemWrite),
    .MemtoReg(MemtoReg),
    .RegWrite(RegWrite),
    .RegDst(RegDst),
    .on_sign(on_sign),
    
    .slt_sign(slt_sign),


    // for j, jr,jal
    .jump(jump),
    .jump_src(jump_src),
    .jump_address(jump_address),
    .jal(jal),

    // for bne,beq
    .branch(branch),
    .branch_address(branch_address),
    .ALUOut(ALUOut),
    .ALUOut_2(ALUOut_Reg[32*2-1:32]),
    .ALUOut_3(ALUOut_Reg[32*3-1:32*2]),


    // for hazard_1
    .hazard_1(hazard_1),

    // for hazard_2
    .hazard_2(hazard_2),
    .RegWrite_2(RegWrite_Reg[0]),
    .RegDst_2(RegDst_Reg[1:0]),
    .rd_2(rd_Reg[4:0]),
    .rt_2(rt_Reg[4:0]),

    // for hazard_3
    .hazard_3(hazard_3),
    .RegWrite_3(RegWrite_Reg[1]),
    .RegDst_3(RegDst_Reg[3:2]),
    .rd_3(rd_Reg[9:5]),
    .rt_3(rt_Reg[9:5]),

    // for MEM-WB hazard
    .Mem_WB_hazard(Mem_WB_hazard),

    // for lw hazard
    .lw_hazard(lw_hazard),
    .MemtoReg_2(MemtoReg_Reg[0]),
    .MemtoReg_3(MemtoReg_Reg[1]),

    // for jr hazard
    .jr_hazard(jr_hazard),

    // the parameter for the WB stage
    .rd_4(rd_Reg[5*2-1:5*1]),
    .rt_4(rt_Reg[5*2-1:5*1]),
    .ALUOut_4(ALUOut_Reg[32*2-1:32]),
    .RegDst_4(RegDst_Reg[2*2-1:2*1]),
    .MemOut_4(MemOut),
    .MemtoReg_4(MemtoReg_Reg[1]),
    .RegWrite_4(RegWrite_Reg[1]),
    .on_sign_4(on_Reg[1]),

    .PCPlus4_4(PCPlus4_Reg[32*3-1:32*2]),
    .jal_4(jal_Reg[1])

);

wire [31:0] ALUOut;
EX ex(
    .clock(clock),
    .rs_data(rs_data),
    .rt_data(rt_data),
    .ALUControl(ALUControl),
    .sa(sa),
    .shift(shift),
    .SignImme(SignImme),
    .ALUSrc(ALUSrc),

    .slt_sign(slt_sign),

    .ALUOut(ALUOut),

    .hazard_1(hazard_1),
    .hazard_2(hazard_2),
    .ALUOut_2(ALUOut_Reg[32*2-1:32]),
    .hazard_3(hazard_3),
    .ALUOut_3(ALUOut_Reg[32*3-1:32*2])


);

wire [31:0] MemOut;

MainMemory MEM(
    .clock(clock),
    .MemWrite(MemWrite_Reg[0]),
    .MemtoReg(MemtoReg_Reg[0]),
    .rt_data(rt_data_Reg[31:0]),
    .ALUOut(ALUOut),

    .MemOut(MemOut),
    .on_sign(on_Reg[1]),

    // for MEM-WB hazard
    .ALUOut_2(ALUOut_Reg[32*2-1:32]),
    .ALUOut_3(ALUOut_Reg[32*3-1:32*2]),
    .ALUOut_4(ALUOut_Reg[32*4-1:32*3]),
    .Mem_WB_hazard(Mem_WB_hazard_Reg[1:0]),

    .stall_times(lw_hazard)

);










    always @(posedge clock) begin

        // pipeline register passing
        // PCPlus4_Reg <= {PC};
        PCF_Reg <= {PCF_Reg[32*3-1:0],PCF};
        rs_data_Reg <= {rs_data_Reg[32*3-1:0],rs_data};
        rt_data_Reg <= {rt_data_Reg[32*3-1:0],rt_data};
        rd_Reg <= {rd_Reg[5*3-1:0],rd};
        rt_Reg <= {rt_Reg[5*3-1:0],rt};
        ALUControl_Reg <= {ALUControl_Reg[4*3-1:0],ALUControl};
        ALUOut_Reg <= {ALUOut_Reg[32*3-1:32*1],ALUOut,{32{1'b0}}};
        sa_Reg <= {sa_Reg[5*3-1:0],sa};;
        shift_Reg <= {shift_Reg[2:0],shift};
        SignImme_Reg <= {SignImme_Reg[32*3-1:0],SignImme};
        ALUSrc_Reg <= {ALUSrc_Reg[2:0],ALUSrc};
        rt_Reg <= {rt_Reg[5*3-1:0],rt };

        if (lw_hazard == 2'b00)
        MemWrite_Reg <= {MemWrite_Reg[2:0],MemWrite};
        else MemWrite_Reg <= {MemWrite_Reg[2:0],1'b0};


        MemtoReg_Reg <= {MemtoReg_Reg[2:0],MemtoReg};



        MemOut_Reg <= {MemOut_Reg[32*3-1:0],MemOut};
        RegDst_Reg <= {RegDst_Reg[2*3-1:0],RegDst};
        on_Reg <= {on_Reg[2:0],on_sign};

        if (lw_hazard == 2'b00)
        RegWrite_Reg <= {RegWrite_Reg[2:0],RegWrite};
        else 
        RegWrite_Reg <= {RegWrite_Reg[2:0],1'b0};
        
        PCPlus4_Reg <= {PCPlus4_Reg[32*3-1:0],PCPlus4};
        jal_Reg <= {jal_Reg[2:0],jal};

        // if (stall_times == 2'b00 )
            stall_times = lw_hazard;
        // if (stall_times != 2'b00)
        // begin
        // stall_times <= stall_times -1'b1;
        // end

        Mem_WB_hazard_Reg <= {Mem_WB_hazard_Reg[2*3-1:0],Mem_WB_hazard};


    end 









endmodule










// the IF stage
module InstructionRAM(clock,instruction,PCPlus4 ,jump,jump_address,jump_src,rs_data, branch,branch_address,stall_times,jr_hazard,ALUOut,ALUOut_2,ALUOut_3);


    input  clock;
    // PCSrcM: 0 : PC => PC + 4, 1 :  PC => PCBranchM

    input jump;
    // Jump: 0: no effect, 1: jump to the given address
    input jump_src;
    // 1 only for jr

    input [31:0] jump_address;
    input [31:0] rs_data;
    input branch;
    input [31:0] branch_address;

    input [1:0] stall_times;

    input [1:0] jr_hazard;
    input [31:0] ALUOut,ALUOut_2,ALUOut_3;


    output reg [31:0] instruction;      
    output reg [31:0] PCPlus4;

    reg [31:0] RAM [0:512-1];
    reg [31:0] PC;

    initial begin
    $readmemb("instuctions.bin",RAM);
    PC = 32'b00000000000000000000000000000000;
    instruction = 32'b00000000000000000000000000000000;

    end





    always @(posedge clock) begin

        // $display("%b",stall_times);
        PCPlus4 <= PC+4;
        if (jump)
        begin
            if (jump_src)
            begin
                if (jr_hazard == 2'b01) 
                    begin
                    if (stall_times == 2'b00)
                        begin
                        instruction <= RAM[ALUOut >>2];
                        PC <= ALUOut+4;
                        end
                    else 
                        begin
                        PC <= PC;
                        instruction <= RAM[(PC-12)>>2];
                        end
                    end
                else if (jr_hazard == 2'b10)
                    begin
                        if (stall_times == 2'b00)
                        begin
                        instruction <= RAM[ALUOut_2 >>2];
                        PC <= ALUOut_2 + 4;
                        end
                    else 
                        begin
                        PC <= PC;
                        instruction <= RAM[(PC-12)>>2];
                        end
                    end
                else if (jr_hazard == 2'b11)
                    begin
                        if (stall_times == 2'b00)
                        begin
                        instruction <= RAM[ALUOut_3 >>2];
                        PC <= ALUOut_3 + 4;
                        end
                    else 
                        begin
                        PC <= PC;
                        instruction <= RAM[(PC-12)>>2];
                        end
                    end
                else 
                    begin
                        if (stall_times == 2'b00)
                        begin
                        instruction <= RAM[rs_data >>2];
                        PC <= rs_data+4;
                        end
                    else 
                        begin
                        PC <= PC;
                        instruction <= RAM[(PC-12)>>2];
                        end
                    end
            end
            else 
                begin
                instruction <= RAM[jump_address >>2];
                if (stall_times == 2'b00)
                PC <= jump_address+4;
                end
        end
        else if (branch) begin

            if (stall_times == 2'b00)
                begin
                instruction <= RAM[(branch_address+PC-8) >> 2];
                PC <= PC+branch_address-4;
                end
            else
                begin 
                    PC <= PC;                instruction <= RAM[(PC-12) >> 2];
                end

        end
        else begin

        if (stall_times == 2'b00)
        begin
        PC <= PC+4;
        instruction <= RAM[PC>>2];
        end
        else 
        begin
        instruction <= RAM[(PC-8)>>2];
        end
        end
        // $display("%b",instruction," ",PC >>2, " stall ",stall_times, " jr hazard ",jr_hazard );
    end

    
endmodule

























// the register part
// two stage is excuted within this module

// first the ID stage
// instruction decodd and register file read
// generate the control signals 

// second the WB stage


module Register(clock,instruction_in, rs_data,rt_data,rd_out,ALUControl,sa_out,shift,SignImme,ALUSrc,MemWrite,MemtoReg,RegWrite,RegDst,rt_out,

on_sign,slt_sign,

jump,jump_address,jump_src,jal,

branch,branch_address, ALUOut,ALUOut_2,ALUOut_3,

hazard_1,

hazard_2,RegWrite_2,RegDst_2,rd_2,rt_2,

hazard_3,RegWrite_3,RegDst_3,rd_3,rt_3,

Mem_WB_hazard,

lw_hazard,MemtoReg_2,MemtoReg_3,

jr_hazard,

rd_4,rt_4,ALUOut_4,RegDst_4,MemOut_4,MemtoReg_4,RegWrite_4,

on_sign_4,

PCPlus4_4,jal_4

);
input clock;
input [31:0] instruction_in;

// for hazard 2
input RegWrite_2;
input [1:0] RegDst_2;
input [4:0] rd_2;
input [4:0] rt_2;

output reg [1:0] hazard_2;

// for hazard 3

input RegWrite_3;
input [1:0] RegDst_3;
input [4:0] rd_3;
input [4:0] rt_3;

output reg [1:0] hazard_3;

output reg [1:0] Mem_WB_hazard;

input MemtoReg_2;
input MemtoReg_3;
output reg [1:0] lw_hazard;

output reg [1:0] jr_hazard;

output reg on_sign;

output reg [3:0] ALUControl; // the 4 bit control signal tells ALU what to do

output reg [31:0]  rs_data,rt_data;

output reg [31:0] SignImme;

output reg [4:0] rd_out,sa_out,rt_out;

output reg ALUSrc; // the control signal tells ALU which data to use
output reg MemtoReg; // only true for lw instuction
output reg MemWrite; // only true for sw instruction


output reg RegWrite;
// the register is written if RegWrite is 1
output reg [1:0]RegDst;
// if RegDst is 00, the writing destination for register is RT
// if RegDst is 01, the writing destination for register is RD
// if RegDst is 11, the writing destination for register is 31 register, this code is unique for JAL instruction

output reg shift;
// if shift is 1, the sa is used instead of rs

output reg slt_sign;
// only on for slt instruction

// the jump is 1 for jump instructions
output reg jump;
output reg [31:0] jump_address;

// only for jr instuction
output reg jump_src;

// only for jal instruction
output reg jal;

// if branch is 1, take the branch
output reg branch;
output reg [31:0] branch_address;
input [31:0] ALUOut,ALUOut_2,ALUOut_3;

output reg [1:0] hazard_1;




// some local variables
reg [5:0] opcode, func_code;
reg [4:0] rs,rt,rd,sa;
reg [15:0] Imme;

// initialize the registers
reg [31:0] Registers [31:0];
integer i;

initial begin
    for (i = 0 ;i < 32; i++)begin
    Registers[i] = 32'b00000000000000000000000000000000;
    end
    lw_hazard = 2'b00;
end 

reg [25:0] target;
reg [31:0] instruction;

reg [31:0] compare_1;
reg [31:0] compare_2;

//————————————————————————————————————————————————————————————————————————————————————————————
// the always of IF stage

always @(posedge clock) begin



if (jump || branch) instruction =32'b00000000000000000000000000000000;
else instruction = instruction_in;

// if instruction is 32'hffffffff, on_sign is 1


if (instruction == 32'hffffffff) on_sign <= 1;
else on_sign <= 0;




// first, take the instruction apart
opcode = instruction[31:26];
func_code = instruction[5:0];

// rs_out <= instruction[25:21];
rt_out <= instruction[20:16];
sa_out <= instruction[10:6];
rd_out <= instruction[15:11];

rs = instruction[25:21];
rt = instruction[20:16];
sa = instruction[10:6];
rd = instruction[15:11];

// do the immediate signed extend
Imme = instruction[15:0];

SignImme <= (Imme[15] == 1'b1)? {{16{1'b1}},Imme} : {{16{1'b0}},Imme};
// $display("Imme",Imme);
// $display("SignImme",SignImme);

// calculate the ALU Control signal
if ((opcode == 6'b001100) || ( opcode == 6'b000000 && func_code == 6'b100100))
ALUControl <= 4'b0000; // and
else if ((opcode == 6'b001101) || ( opcode == 6'b000000 && func_code == 6'b100101))
ALUControl <= 4'b0001; // OR
else if ((opcode == 6'b100011)||(opcode == 6'b101011)||(opcode == 6'b001000)||(opcode == 6'b001001) || ( opcode == 6'b000000 && func_code == 6'b100000)||( opcode == 6'b000000 && func_code == 6'b100001))
ALUControl <= 4'b0010; // add

else if (( opcode == 6'b000000 && func_code == 6'b100010)||( opcode == 6'b000000 && func_code == 6'b100011))
ALUControl <= 4'b0011; // sub
else if (( opcode == 6'b000000 && func_code == 6'b101010))
ALUControl <= 4'b0100; // set on less than

else if (func_code == 6'b100111 && opcode == 6'b000000)
ALUControl <= 4'b0101; // nor

else if (opcode == 6'b001110 ||( opcode == 6'b000000 && func_code == 6'b100110))
ALUControl <= 4'b0110; // xor

else if (( opcode == 6'b000000 && func_code == 6'b000000)||( opcode == 6'b000000 && func_code == 6'b000100))
ALUControl <= 4'b0111; // shift left

else if (( opcode == 6'b000000 && func_code == 6'b000010)||( opcode == 6'b000000 && func_code == 6'b000110))
ALUControl <= 4'b1000; // shift right logically

else if (( opcode == 6'b000000 && func_code == 6'b000011)||( opcode == 6'b000000 && func_code == 6'b000111))

ALUControl <= 4'b1001; // shift right arithmatically

else ALUControl <= 4'b1111; // if the instruction do not use the ALU
// that's all of the ALU control signal


// generate the ALUSrc Signal 
// if the signal is 01, the immediate take RT's place in ALU   
if (opcode == 6'b100011 // lw
|| opcode == 6'b101011 // sw
|| opcode == 6'b001000 // addi
|| opcode == 6'b001001 // addiu
|| opcode == 6'b001100 // andi
|| opcode == 6'b001101 // ori
|| opcode == 6'b001110 // xori
)
ALUSrc <= 1'b1;
else ALUSrc <= 1'b0;



// generate the MemtoReg signal
if (opcode == 6'b100011) 
begin
MemtoReg <= 1'b1; //lw

end
else MemtoReg <= 1'b0;

// generate the MemWrite signal
if (opcode == 6'b101011 ) MemWrite <= 1'b1; //sw
else MemWrite <= 1'b0;



// generate RegWrite
if ( 
// jr
(opcode == 6'b000000 & func_code == 6'b001000)
// beq
|| opcode == 6'b000100  
// bne
|| opcode == 6'b000101 
// j
|| opcode == 6'b000010
//sw
|| opcode == 6'b101011
// initial stage
|| instruction == 32'b00000000000000000000000000000000
|| instruction == 32'b11111111111111111111111111111111
)
 RegWrite <= 1'b0;
else  RegWrite <= 1'b1; 

// $display("%b",opcode,ALUSrc,RegWrite);

// $display(opcode,"   ",RegWrite);

// generate RegDst
if (opcode == 6'b000000) // R type
RegDst <= 2'b01;
else if (opcode == 6'b000011)// jal
RegDst <= 2'b11;
else RegDst <= 2'b00; // I type

// generate shift
if ((opcode == 6'b000000 && func_code == 6'b000000)|| (opcode == 6'b000000 && func_code == 6'b000011)|| (opcode == 6'b000000 && func_code == 6'b000010))
shift <= 1'b1;
else shift <= 1'b0;

// generate slt_sign
if (opcode == 6'b000000 && func_code == 6'b101010)
slt_sign <= 1'b1;
else slt_sign <= 1'b0;

//generate jump and jump address
if ((opcode == 6'b000000 && func_code == 6'b001000) // jr
|| opcode == 6'b000010 // j
|| opcode == 6'b000011 // jal
)
jump <= 1'b1;
else jump <= 1'b0;



target = instruction[25:0];
jump_address <= {4'b0000,target,2'b00};


// generate jump src
if ((opcode == 6'b000000 && func_code == 6'b001000)) // jr
jump_src <=1;
else jump_src <= 0;

// generate jal
if (opcode == 6'b000011)
jal <=1;
else jal <= 0;








// generate brach with hazard

// first, get comapre_1,compare_2



if (opcode == 6'b000100 && Registers[rs] == Registers[rt])
branch <= 1'b1;
else if (opcode == 6'b000101 && Registers[rs] != Registers[rt])
branch <= 1'b1;
else branch <= 1'b0;

branch_address <= (Imme[15] == 1'b1) ? {{14{1'b1}},Imme,{2{1'b0}}} : {{14{1'b0}},Imme,{2{1'b0}}} ;


// $display("%b",opcode,branch,comapre_1,comapre_2);





// then, fetch data from registers
rs_data <= Registers[rs];
rt_data <= Registers[rt];

// $display("%b",instruction);

// detect hazard_1;



if (RegWrite) begin
    if (RegDst == 2'b01 && (rd_out == rs))
        hazard_1 <= 2'b01; // rs will be forward
    else if (RegDst == 2'b01 && (rd_out == rt) && opcode == 6'b000000 )
        hazard_1 <= 2'b11; // rt will be forward
    else if (RegDst ==2'b00 && (rt_out == rs))
        hazard_1 <= 2'b01;
    else if (RegDst == 2'b00 && (rt_out == rt) && opcode == 6'b000000)
        hazard_1 <= 2'b11;
    else hazard_1 <= 2'b00;
end
else 
hazard_1 <= 2'b00;

// detect hazard_2;
if (RegWrite_2) begin
    if (RegDst_2 == 2'b01 && (rd_2 == rs))
        hazard_2 <= 2'b01; // rs will be forward
    else if (RegDst_2 == 2'b01 && (rd_2 == rt) && opcode == 6'b000000 )
        hazard_2 <= 2'b11; // rt will be forward
    else if (RegDst_2 ==2'b00 && (rt_2 == rs))
        hazard_2 <= 2'b01;
    else if (RegDst_2 == 2'b00 && (rt_2 == rt) && opcode == 6'b000000)
        hazard_2 <= 2'b11;
    else hazard_2 <= 2'b00;
end
else 
hazard_2 <= 2'b00;

// detect hazard_3

if (RegWrite_3) begin
    if (RegDst_3 == 2'b01 && (rd_3 == rs))
        hazard_3 <= 2'b01; // rs will be forward
    else if (RegDst_3 == 2'b01 && (rd_3 == rt) && opcode == 6'b000000 )
        hazard_3 <= 2'b11; // rt will be forward
    else if (RegDst_3 ==2'b00 && (rt_3 == rs))
        hazard_3 <= 2'b01;
    else if (RegDst_3 == 2'b00 && (rt_3 == rt) && opcode == 6'b000000)
        hazard_3 <= 2'b11;
    else hazard_3 <= 2'b00;
end
else 
hazard_3 <= 2'b00;

// $display(hazard_1,hazard_2,hazard_3);

// detect MEM-WB hazard

if (opcode == 6'b101011 && RegWrite == 1'b1 && RegDst == 2'b01 && rt == rd_out)
Mem_WB_hazard <= 2'b01;
else if (opcode == 6'b101011 && RegWrite == 1'b1 && RegDst == 2'b00 && rt == rt_out)
Mem_WB_hazard <= 2'b01;
else if (opcode == 6'b101011 && RegWrite_2 == 1'b1 && RegDst_2 == 2'b01 && rt == rd_2)
Mem_WB_hazard <= 2'b10;
else if (opcode == 6'b101011 && RegWrite_2 == 1'b1 && RegDst_2 == 2'b00 && rt == rt_2)
Mem_WB_hazard <= 2'b10;
else if (opcode == 6'b101011 && RegWrite_3 == 1'b1 && RegDst_3 == 2'b01 && rt == rd_3)
Mem_WB_hazard <= 2'b11;
else if (opcode == 6'b101011 && RegWrite_3 == 1'b1 && RegDst_3 == 2'b00 && rt == rt_3)
Mem_WB_hazard <= 2'b11;
else Mem_WB_hazard <= 2'b00;

// $display(Mem_WB_hazard);

// detect lw hazard
if (lw_hazard ==2'b00)
begin
if (MemtoReg && (rs == rt_out || rt == rt_out) && opcode == 6'b000000 )
begin
    lw_hazard <= 2'b11;
end
else if (MemtoReg && (rs == rt_out) && opcode != 6'b000000 )
    lw_hazard <= 2'b11;
else if (MemtoReg_2 && (rs == rt_2 || rt == rt_2) && opcode == 6'b000000 )
    lw_hazard <= 2'b10;
else if (MemtoReg_2 && (rs == rt_2) && opcode != 6'b000000 )
    lw_hazard <= 2'b10;
else if (MemtoReg_3 && (rs == rt_3 || rt == rt_3) && opcode == 6'b000000 )
    lw_hazard <= 2'b01;
else if (MemtoReg_3 && (rs == rt_3) && opcode != 6'b000000 )
    lw_hazard <= 2'b01;
// lw and sw conflicts
else if (opcode == 6'b101011 && rt == rt_out && MemtoReg)
    lw_hazard <= 2'b11;
else if (opcode == 6'b101011 && rt == rt_2 && MemtoReg_2)
    lw_hazard <= 2'b10;
else if (opcode == 6'b101011 && rt == rt_3&& MemtoReg_3)
    lw_hazard <= 2'b01;
else lw_hazard <= 2'b00;
end
else 
lw_hazard <= lw_hazard -1;


// generate jr hazard
if (lw_hazard == 2'b00)
begin
if (opcode == 6'b000000 && func_code == 6'b001000)
begin
    if (RegWrite && RegDst == 2'b00 && rt_out == rs)
        jr_hazard <= 2'b01;
    else if (RegWrite && RegDst == 2'b01 && rd_out == rs)
        jr_hazard <= 2'b01;
    else if (RegWrite_2 && RegDst_2 == 2'b00 && rt_2 == rs)
    jr_hazard <= 2'b10;
else if (RegWrite_2 && RegDst_2 == 2'b01 && rd_2 == rs)
    jr_hazard <= 2'b10;
    else if (RegWrite_3 && RegDst_3 == 2'b00 && rt_3 == rs)
    jr_hazard <= 2'b11;
else if (RegWrite_3 && RegDst_3 == 2'b01 && rd_3 == rs)
    jr_hazard <= 2'b11;
    else jr_hazard <= 2'b00;
end
else 
jr_hazard <= 2'b00;
end
else jr_hazard <= 2'b00;





// $display("%b",instruction," hi ","%b",MemtoReg_2," rt ",rt," rt out 2 ", rt_2,"  lw hazard  ",lw_hazard, "opcode %b",opcode);

end 


//————————————————————————————————————————————————————————————————————————————————————————————
// the always for WB stage


input [31:0] MemOut_4;
input [31:0] ALUOut_4;
input [4:0] rt_4,rd_4;
input [31:0] PCPlus4_4;
input jal_4;

input MemtoReg_4;
input [1:0] RegDst_4;
input RegWrite_4;

input on_sign_4;

reg [31:0] content;
reg [4:0] WritDst;

output reg on;

always @(posedge clock) begin

    


    // first, decide the writing content
    if (jal_4) 
    begin
    content = PCPlus4_4;
    end
    else 
    begin
        if (MemtoReg_4) content = MemOut_4;
        else content = ALUOut_4;
    end


    // $display(RegWrite_4);

    // then, decide the writing destination
    if (RegDst_4 == 2'b00)
    begin
        WritDst = rt_4;
    end
    else if (RegDst_4 == 2'b01)
        WritDst = rd_4;
    else if (RegDst_4 == 2'b11)
    begin
        WritDst = 31;
    end

    // write to the register
    if (RegWrite_4)
    begin
        // $display(WritDst,$signed(content)," ",RegDst_4);
        Registers[WritDst] <= content;
    end

end



endmodule

















module EX(clock,rs_data,rt_data,sa,SignImme,shift,ALUSrc,ALUOut,ALUControl, slt_sign, hazard_1,hazard_2,ALUOut_2,hazard_3,ALUOut_3);
input [31:0] rs_data,rt_data;
input [31:0] SignImme;
input [3:0] ALUControl;
input [4:0] sa;
input shift;
input ALUSrc;

input slt_sign;

input clock;
// input [31:0] PCPlus4;

input [1:0] hazard_1;
input [1:0] hazard_2;
input [31:0] ALUOut_2;

input [1:0] hazard_3;
input [31:0] ALUOut_3;

output reg [31:0]  ALUOut;
// output reg [31:0] PCBranch;


reg [31:0] input_1;
reg [31:0] input_2;
reg [31:0] result;


always @(posedge clock) begin

// $display("IN EX",SignImme);

// first, get the correct input
if (shift == 1) input_1 = sa;
else if (hazard_1 == 2'b01) 
begin
input_1 = ALUOut;
// $display(ALUOut);
end
else if (hazard_1 !=2'b01 && hazard_2 == 2'b01) input_1 = ALUOut_2;
else if (hazard_1 !=2'b01 && hazard_2 != 2'b01 && hazard_3 == 2'b01)
begin 
input_1 = ALUOut_3;
// $display($signed(ALUOut_3));
end
else input_1 = rs_data;  

if (ALUSrc == 1) input_2 = SignImme;
else if (hazard_1 == 2'b11) input_2 = ALUOut;
else if (hazard_1 != 2'b11 && hazard_2 == 2'b11 ) 
begin
input_2 = ALUOut_2;
// $display(ALUOut_2);
end
else if (hazard_1 != 2'b11 && hazard_2 != 2'b11 && hazard_3 == 2'b11) input_1 = ALUOut_3;
else input_2 = rt_data;

// $display($signed(input_1),$signed(input_2));
// do the ALU calculation
case(ALUControl)
    4'b0000: // and
    begin
    if (ALUSrc) ALUOut <= input_1 & input_2[15:0];
    else ALUOut <= input_1 & input_2;
    end
    4'b0001: // or
    begin
    if (ALUSrc) ALUOut <= input_1 | input_2[15:0];
    else ALUOut <= input_1 | input_2;
    end
    4'b0010: // add
    ALUOut <= input_1 + input_2;
    4'b0011: // sub
    ALUOut <= input_1 - input_2;
    4'b0100:// set on less than
    begin
    result = input_1 - input_2;
    ALUOut <= {{31{1'b0}},{result[31]}} ;
    end
    4'b0101: // nor
    begin
    ALUOut <= ~(input_1 | input_2);
    end
    4'b0110: // xor
    begin
    if (ALUSrc) ALUOut <= input_1 ^ input_2[15:0];
    else ALUOut <= input_1 ^ input_2;
    end
    4'b0111: // shift left
    ALUOut <= input_2 << input_1[4:0];
    4'b1000: // shift right logically
    ALUOut <= input_2 >> input_1[4:0];
    4'b1001: // shift right arithmatically
    ALUOut <= $signed(input_2) >>> input_1[4:0];
endcase




// // calculate the branch address
// PCBranch = (SignImme << 2) + PCPlus4;
// $display($signed(ALUOut));

end // the end of always
endmodule






















// the MEM stage
module MainMemory (clock,MemWrite,MemtoReg,rt_data,ALUOut,MemOut,on_sign,ALUOut_2,Mem_WB_hazard,stall_times,ALUOut_3,ALUOut_4);
input clock; 

input [31:0] rt_data;
input [31:0] ALUOut;

input MemWrite;
input MemtoReg;
input on_sign;

input [31:0] ALUOut_4,ALUOut_2,ALUOut_3;
input [1:0] Mem_WB_hazard;
input [1:0] stall_times;


output reg [31:0] MemOut;

// initialize the memory
reg [31:0] DATA_RAM [0:512-1];

reg [16383:0] ram_init;
integer i;
integer f;

    initial begin
    // this kind of initialization may be slower?
        ram_init = {16384{1'b0}};
    for (i=0; i < 512; i = i + 1) begin
        DATA_RAM[512-1-i] = ram_init[i*32+:32];
    end
    f = $fopen("result_5.txt","w");

    end


    always @(posedge clock) begin
        
        // $display(Mem_WB_hazard,MemWrite);
        // $display($signed(ALUOut >> 2),$signed(ALUOut_2)); 
        if (MemWrite) begin //sw;
            if (Mem_WB_hazard == 2'b01)
            begin
            DATA_RAM[ALUOut >> 2] = ALUOut_2;
            MemOut <= {{32{1'b0}}};
            end
            else if (Mem_WB_hazard == 2'b10)
            begin
            DATA_RAM[ALUOut >> 2] = ALUOut_3;
            MemOut <= {{32{1'b0}}};
            end
            else if (Mem_WB_hazard == 2'b11)
            begin
            DATA_RAM[ALUOut >> 2] = ALUOut_4;
            MemOut <= {{32{1'b0}}};
            end
            else begin
            DATA_RAM[ALUOut >> 2] = rt_data;
            // $display("sw",ALUOut >>2, rt_data);
            MemOut <= {{32{1'b0}}};
            end

        end
        else if(MemtoReg) // lw
        begin
            MemOut <= DATA_RAM[ALUOut >> 2];
        end
        else MemOut <= {32{1'b0}};




        // $display(on_sign);
        if (on_sign && stall_times == 2'b00)
        begin
            for (i=0; i < 512; i = i + 1) begin
            $fwrite(f,"%b\n",DATA_RAM[i]);
            end
            $fclose(f);
            $finish;
        end
    end

endmodule


